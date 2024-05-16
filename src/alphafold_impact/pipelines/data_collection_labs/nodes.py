"""This module contains the nodes used in the data collection pipeline for labs.

The nodes are used to collect papers from the OpenAIRE API, extract candidate
authors from the papers, and fetch publications for the candidate authors. 
"""

import logging
import time
from typing import Dict, List, Tuple, Generator, Union, Iterator
import requests
from requests.adapters import HTTPAdapter, Retry
from kedro.io import AbstractDataset
import pandas as pd
from joblib import Parallel, delayed
from ..data_collection_oa.nodes import collect_papers  # pylint: disable=E0402

logger = logging.getLogger(__name__)


def get_candidate_authors(
    data: pd.DataFrame, baseline: bool = False
) -> List[Tuple[str, str]]:
    """
    Get the last authors of the given papers.

    Args:
        data (pd.DataFrame): The input DataFrame containing the papers.

    Returns:
        List[Tuple[str, str]]: The list of last authors.
    """
    # threshold = 10 if baseline else 1

    # get the author whose third element is the last author
    author_data = (
        data.drop_duplicates(subset=["id"])
        .explode("authorships")
        .assign(
            author=lambda x: x["authorships"].apply(
                lambda y: y[0] if y is not None else None
            ),
            institution=lambda x: x["authorships"].apply(
                lambda y: y[1] if y is not None else None
            ),
            position=lambda x: x["authorships"].apply(
                lambda y: y[2] if y is not None else None
            ),
        )
        .dropna(subset=["author"])
        .drop_duplicates(subset=["id", "author", "institution", "position"])
        .reset_index(drop=True)
    )

    # drop A9999999999
    author_data = author_data[~(author_data["author"] == "A9999999999")]

    # drop empty string institution
    author_data = author_data[~(author_data["institution"] == "")]

    if baseline:
        # filter to keep only relevant authors (5> as last, or first author)
        author_institution_counts = author_data.groupby(
            ["author", "institution", "position"]
        )["id"].nunique()
        filtered_pairs = author_institution_counts[
            (author_institution_counts >= 5)
            & (
                author_institution_counts.index.get_level_values("position").isin(
                    ["first", "last"]
                )
            )
        ]
        baseline_authors = [
            (author, institution)
            for author, institution, _ in list(filtered_pairs.index)
        ]

        return list(set(baseline_authors))

    # instead return all author, institution pairs when position is last
    alphafold_authors = (
        author_data[author_data["position"] == "last"]
        .groupby(["author", "institution"])["id"]
        .nunique()
    )

    return list(alphafold_authors.index)


def merge_candidate_authors(data: AbstractDataset) -> List[Tuple[str, str]]:
    """
    Merges the candidate authors from the given dataset.

    Args:
        data (AbstractDataset[Dict[str, List[str]]]): The dataset containing candidate authors.

    Returns:
        List[Tuple[str, str]]: A list of tuples representing the merged candidate authors.

    """
    candidates = [item for sublist in data.values() for item in sublist()]

    return list(set(tuple(i) for i in candidates))


def fetch_author_publications(
    author_ids: List[str],
    from_publication_date: str,
    api_config: Dict[str, str],
) -> Generator[Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]], None, None]:
    """
    Fetches publications for a list of authors based on their IDs and publication date.

    Args:
        author_ids (List[str]): List of author IDs.
        from_publication_date (str): Starting publication date in the format 'YYYY-MM-DD'.
        api_config (Dict[str, str]): API configuration dictionary containing 'mailto' and
            'perpage' values.

    Yields:
        Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]]: A generator that yields a
            tuple of dictionaries. The first dictionary contains dataframes with publication
            information, where the keys are slice IDs. The second dictionary contains lists
            of publication dictionaries, where the keys are slice IDs.

    """
    author_ids = [author[0] for author in author_ids]

    # create batches of 50 authors
    author_batches = [
        "|".join(author_ids[i : i + 25]) for i in range(0, len(author_ids), 25)
    ]

    filter_batches = [[from_publication_date, batch] for batch in author_batches]

    # slice to create parallel jobs that produce slices
    slices = [filter_batches[i : i + 250] for i in range(0, len(filter_batches), 250)]

    logger.info("Fetching papers for %d author batches", len(slices))

    for i, slice_ in enumerate(slices):
        logger.info("Processing batch number %d / %d", i + 1, len(slices))

        slice_results = Parallel(n_jobs=8, backend="loky", verbose=10)(
            delayed(collect_papers)(
                oa_ids=batch,
                mailto=api_config["mailto"],
                perpage=api_config["perpage"],
                filter_criteria=["from_publication_date", "author.id"],
                slice_keys=True,
                eager_loading=True,
                skip_preprocess=True,
                sample_size=200 * len(batch[1].split("|")),
            )
            for batch in slice_
        )

        slice_papers = [
            paper
            for sublist in [v for d in slice_results for k, v in d.items()]
            for paper in sublist
        ]

        # for each listed dict, keep "id", "authorships", "publication_date", "cited_by_count"
        slice_papers = [
            {
                "id": paper["id"].replace("https://openalex.org/", ""),
                "authorships": paper["authorships"],
                "publication_date": paper["publication_date"],
                "cited_by_count": paper["cited_by_count"],
            }
            for paper in slice_papers
        ]

        yield {f"s{i}": slice_papers}


def _works_generator(
    mailto: str,
    perpage: str,
    oa_id: Union[str, List[str]],
    session: requests.Session,
) -> Iterator[list]:
    """Creates a generator that yields a list of works from the OpenAlex API based on a
    given work ID.

    Args:
        mailto (str): The email address to use for the API.
        perpage (str): The number of results to return per page.
        oa_id (Union[str, List[str]): The work ID to use for the API.
        session (requests.Session): The requests session to use.

    Yields:
        Iterator[list]: A generator that yields a list of works from the OpenAlex API
        based on a given work ID.
    """
    cursor = "*"

    cursor_url = (
        f"https://api.openalex.org/authors?search={oa_id}"
        f"&filter=affiliations.institution.country_code:US"
        f"&mailto={mailto}&per-page={perpage}&cursor={{}}"
    )

    try:
        # make a call to estimate total number of results
        response = session.get(cursor_url.format(cursor), timeout=20)
        data = response.json()

        while response.status_code == 429:  # needs testing (try with 200)
            logger.info("Waiting for 1 hour...")
            time.sleep(3600)
            response = session.get(cursor_url.format(cursor), timeout=20)
            data = response.json()

        logger.info("Fetching data for %s", oa_id[:50])
        total_results = data["meta"]["count"]
        num_calls = total_results // int(perpage) + 1
        logger.info("Total results: %s, in %s calls", total_results, num_calls)
        while cursor:
            response = session.get(cursor_url.format(cursor), timeout=20)
            data = response.json()
            results = data.get("results")
            cursor = data["meta"].get("next_cursor", False)
            yield results

    except Exception as e:  # pylint: disable=broad-except
        logger.error("Error fetching data for %s: %s", oa_id, e)
        yield []


def fetch_ids_for_author(
    oa_id: Union[str, List[str]], mailto: str, perpage: str, **kwargs
) -> List[dict]:
    """Fetches all papers cited by a specific work ID."""
    author_ids = []
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=0.3)
    session.mount("https://", HTTPAdapter(max_retries=retries))
    for page, authors in enumerate(
        _works_generator(mailto, perpage, oa_id, session, **kwargs)
    ):
        author_ids.extend(authors)
        logger.info(
            "Fetching page %s. Total papers collected: %s",
            page,
            len(author_ids),
        )

    return author_ids


def get_pi_id(data: pd.DataFrame) -> pd.DataFrame:
    """
    This function fetches the PI IDs from the given data. It extracts the PI name from the
    'principal_investigator' column, fetches the IDs for the PI name, and then extracts the
    PI ID, display name, last known institution display name, last known institution type,
    and last known institution name.

    Args:
        data (pd.DataFrame): The input DataFrame containing the data.

    Returns:
        List[str]: The list of PI IDs.
    """

    # extract the PI name, by splitting and taking left before Ph.D or M.D
    data["pi_name"] = data["principal_investigator"].str.split("Ph.D.|M.D.").str[0]

    # split on ",", and invert order
    data["pi_name"] = data["pi_name"].str.split(",").str[::-1].str.join(" ").str.strip()

    # for each name, fetch ids
    data["pi_oa"] = data["pi_name"].apply(
        lambda x: fetch_ids_for_author(
            x,
            mailto="david.ampudia@nesta.org.uk",
            perpage="100",
        )
    )

    # explode data
    data = data.explode("pi_oa")

    # drop na
    data = data.dropna(subset=["pi_oa"])

    # create relevant columns
    data["pi_id"] = data["pi_oa"].apply(lambda x: x.get("id"))
    data["pi_display_name"] = data["pi_oa"].apply(lambda x: x.get("display_name"))
    data["pi_last_known_institution"] = data["pi_oa"].apply(
        lambda x: x.get("last_known_institution")
    )

    # drop None of pi_last_known_institution
    data = data.dropna(subset=["pi_last_known_institution"])

    # extract name and type
    data["pi_last_known_institution_name"] = data["pi_last_known_institution"].apply(
        lambda x: x.get("display_name")
    )
    data["pi_last_known_institution_type"] = data["pi_last_known_institution"].apply(
        lambda x: x.get("type")
    )

    # keep rows with National in past_last_known_institution_type
    data = data[data["pi_last_known_institution_name"].str.contains("National")]

    # filter for those that have pi_last_known_institution_type government
    data = data[
        data["pi_name"].duplicated(keep=False)
        & data["pi_last_known_institution_type"].str.contains("government")
    ]

    # if still duplicates, take the first
    data = data.drop_duplicates(subset=["pi_name"], keep="first")

    return data


def get_publications_from_labs(
    data: pd.DataFrame,
    from_publication_date: str,
    api_config: Dict[str, str],
) -> Generator[Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]], None, None]:
    """
    Fetches publications for a given list of authors.

    Args:
        data (pd.DataFrame): The input data containing author information.
        from_publication_date (str): The starting publication date to filter the publications.
        api_config (Dict[str, str]): Configuration settings for the API.

    Yields:
        Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]]: A generator that yields a tuple 
            of dictionaries. The first dictionary contains slice keys as keys and a list of papers 
            as values. The second dictionary contains slice keys as keys and a list of authors as 
            values.

    Returns:
        None: This function does not return anything directly. It yields the results using a 
            generator.
    """

    # get unique authors from data
    author_ids = data["author"].unique().tolist()

    logger.info("Fetching publications for %d authors", len(author_ids))

    # # create batches of 50 authors
    # author_batches = [
    #     "|".join(author_ids[i : i + 50]) for i in range(0, len(author_ids), 50)
    # ]

    filter_batches = [[from_publication_date, author_id] for author_id in author_ids]

    # slice to create parallel jobs that produce slices
    slices = [filter_batches[i : i + 250] for i in range(0, len(filter_batches), 250)]

    logger.info("Fetching papers for %d author batches", len(slices))

    for i, slice_ in enumerate(slices):

        logger.info("Processing batch number %d / %d", i + 1, len(slices))

        slice_results = Parallel(n_jobs=8, backend="loky", verbose=10)(
            delayed(collect_papers)(
                oa_ids=batch,
                mailto=api_config["mailto"],
                perpage=api_config["perpage"],
                filter_criteria=["from_publication_date", "author.id"],
                slice_keys=True,
                eager_loading=True,
                skip_preprocess=True,
            )
            for batch in slice_
        )

        slice_papers = {
            author_id: papers
            for (_, author_id), slice_result in zip(slice_, slice_results)
            for papers in slice_result.values()
        }

        logger.info("Processing batch number %d / %d", i + 1, len(slices))
        output = []

        for author_id, children_list in slice_papers.items():

            json_data = [
                {
                    k: v
                    for k, v in item.items()
                    if k
                    in [
                        "id",
                        "ids",
                        "doi",
                        "display_name",
                        "publication_date",
                        "mesh_terms",
                        "cited_by_count",
                        "counts_by_year",
                        "authorships",
                        "concepts",
                        "topics",
                    ]
                }
                for item in children_list
            ]

            # Check if json_data is not empty
            if len(json_data) > 0:
                # Transform to DataFrame and add parent_id
                df = pd.DataFrame(json_data)
                df["pi_id"] = author_id

                # Extract pmid from ids
                df["pmid"] = df["ids"].apply(
                    lambda x: (
                        x.get("pmid").replace("https://pubmed.ncbi.nlm.nih.gov/", "")
                        if x and x.get("pmid")
                        else None
                    )
                )

                # transform to dataframe and add parent_id
                df = pd.DataFrame(json_data)
                df["pi_id"] = author_id

                # extract pmid from ids
                df["pmid"] = df["ids"].apply(
                    lambda x: (
                        x.get("pmid").replace("https://pubmed.ncbi.nlm.nih.gov/", "")
                        if x and x.get("pmid")
                        else None
                    )
                )

                # keep only a list of tuples with "descriptor_ui" and "descriptor_name" for each
                df["mesh_terms"] = df["mesh_terms"].apply(
                    lambda x: (
                        [(c["descriptor_ui"], c["descriptor_name"]) for c in x]
                        if x
                        else None
                    )
                )

                # break atuhorship nested dictionary jsons, create triplets of authorship
                df["authorships"] = df["authorships"].apply(
                    lambda x: (
                        [
                            (
                                (
                                    author["author"]["id"].replace(
                                        "https://openalex.org/", ""
                                    ),
                                    inst["id"].replace("https://openalex.org/", ""),
                                    author["author_position"],
                                )
                                if author["institutions"]
                                else [
                                    author["author"]["id"].replace(
                                        "https://openalex.org/", ""
                                    ),
                                    "",
                                    author["author_position"],
                                ]
                            )
                            for author in x
                            for inst in author["institutions"] or [{}]
                        ]
                        if x
                        else None
                    )
                )

                # create a list of topics
                df["topics"] = df["topics"].apply(
                    lambda x: (
                        [
                            (
                                topic["id"].replace("https://openalex.org/", ""),
                                topic["display_name"],
                                topic["subfield"]["id"].replace(
                                    "https://openalex.org/", ""
                                ),
                                topic["subfield"]["display_name"],
                                topic["field"]["id"].replace(
                                    "https://openalex.org/", ""
                                ),
                                topic["field"]["display_name"],
                                topic["domain"]["id"].replace(
                                    "https://openalex.org/", ""
                                ),
                                topic["domain"]["display_name"],
                            )
                            for topic in x
                        ]
                        if x
                        else None
                    )
                )

                # extract concepts
                df["concepts"] = df["concepts"].apply(
                    lambda x: (
                        [
                            (
                                concept["id"].replace("https://openalex.org/", ""),
                                concept["display_name"],
                            )
                            for concept in x
                        ]
                        if x
                        else None
                    )
                )

                # change doi to remove the url
                df["doi"] = df["doi"].str.replace("https://doi.org/", "")

                # change id to remove the url
                df["id"] = df["id"].str.replace("https://openalex.org/", "")

                # append to output
                output.append(df)

        df = pd.concat(output)

        yield {f"s{i}": df}
