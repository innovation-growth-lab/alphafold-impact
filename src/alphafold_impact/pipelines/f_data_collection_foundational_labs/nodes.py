"""
This is a boilerplate pipeline 'data_collection_sb_labs'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict, List, Tuple, Generator
from kedro.io import AbstractDataset
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from sklearn.preprocessing import MinMaxScaler

from ..a_data_collection_oa.nodes import collect_papers  # pylint: disable=E0402
from .utils import (
    appears_in_three_consecutive_years,
    compute_avg_citation_count,
    compute_sample_publication_count,
    compute_apf,
    explode_author_data,
    fetch_institution_data,
    get_sb_candidate_authors,
    is_likely_pi,
    separate_ct_from_seed,
)

logger = logging.getLogger(__name__)


def get_candidate_authors(
    alphafold_data: pd.DataFrame,
    baseline_data: pd.DataFrame,
    seed_data: pd.DataFrame,
    ct_data: pd.DataFrame,
):
    """
    Get a list of candidate authors from different data sources.

    Args:
        alphafold_data (pd.DataFrame): Dataframe containing AlphaFold data.
        baseline_data (pd.DataFrame): Dataframe containing baseline data.
        seed_data (pd.DataFrame): Dataframe containing seed data.
        ct_data (pd.DataFrame): Dataframe containing CT data.

    Returns:
        list: A list of candidate authors from different data sources.
    """
    # get the last authors from the AlphaFold data
    alphafold_authors = get_sb_candidate_authors(alphafold_data)

    # separate the CT data from the seed data
    ct_long_ai_data, ct_long_pp_data, ct_long_sb_data, other_data = (
        separate_ct_from_seed(baseline_data, seed_data, ct_data)
    )

    # get the last authors from the CT data
    ct_ai_authors = get_sb_candidate_authors(ct_long_ai_data)
    ct_pp_authors = get_sb_candidate_authors(ct_long_pp_data)
    ct_sb_authors = get_sb_candidate_authors(ct_long_sb_data)

    # add the CT authors
    ct_authors = ct_ai_authors + ct_pp_authors + ct_sb_authors

    # get the last authors from the other SB data
    other_authors = get_sb_candidate_authors(other_data)

    # create a unique list of authors
    authors = list(set(alphafold_authors + ct_authors + other_authors))

    return (
        authors,
        alphafold_authors,
        ct_ai_authors,
        ct_pp_authors,
        ct_sb_authors,
        ct_authors,
        other_authors,
    )


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
        api_config (Dict[str, str]): API configuration dictionary containing
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
    slices = [filter_batches[i : i + 25] for i in range(0, len(filter_batches), 25)]

    logger.info("Fetching papers for %d author batches", len(slices))

    for i, slice_ in enumerate(slices):
        if i < 718:
            continue
        logger.info("Processing batch number %d / %d", i + 1, len(slices))

        slice_results = Parallel(n_jobs=8, backend="loky", verbose=10)(
            delayed(collect_papers)(
                oa_ids=batch,
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


def calculate_lab_determinants(
    dict_loader: AbstractDataset,
    alphafold_authors: List,
    ct_ai_authors: List,
    ct_pp_authors: List,
    ct_sb_authors: List,
    other_authors: List,
) -> Generator[Dict[str, pd.DataFrame], None, None]:
    """
    This node calculates the lab determinants for each lab based on the given data.

    Parameters:
        dict_loader (AbstractDataset): The input data loader containing the lab data.
        alphafold_authors (List): The list of AlphaFold authors.
        ct_ai_authors (List): The list of CT AI authors.
        ct_pp_authors (List): The list of CT PP authors.
        ct_sb_authors (List): The list of CT SB authors.
        other_authors (List): The list of other authors.

    Yields:
        Generator[Dict[str, pd.DataFrame], None, None]: A dictionary containing the lab
            determinants for each lab.

    """
    for i, (key, loader) in enumerate(dict_loader.items()):
        logger.info(
            "Processing data from %s. Number %s out of %s", key, i + 1, len(dict_loader)
        )

        article_list = loader()

        json_data = [
            {
                k: v
                for k, v in item.items()
                if k
                in [
                    "id",
                    "doi",
                    "publication_date",
                    "cited_by_count",
                    "authorships",
                ]
            }
            for item in article_list
        ]

        # transform to dataframe and add parent_id
        data = pd.DataFrame(json_data)

        # break authorship nested dictionary jsons, create pipe-delimited string of authorship triplets
        data["authorships"] = data["authorships"].apply(
            lambda x: (
                "|".join([
                    f"{author['author']['id'].replace('https://openalex.org/', '')},{inst['id'].replace('https://openalex.org/', '') if inst else ''},{author['author_position']}"
                    for author in x
                    for inst in author["institutions"] or [{}]
                ]) if x else None
            )
        )

        # for each dictionary, only preserve
        author_data = explode_author_data(data)

        # get mapping between candidates and chain seeds
        candidate_map = create_candidates_map(
            alphafold_authors,
            ct_ai_authors,
            ct_pp_authors,
            ct_sb_authors,
            other_authors,
        )

        # keep if author, institution pair is in candidate_map author institution columns
        author_data = author_data.merge(
            candidate_map[["author", "institution", "seed"]],
            on=["author", "institution"],
            how="inner",
        )

        # lets drop authors whose row count is less than 10 (reduces chances of reassessment)
        author_data = author_data.groupby(["author", "institution"]).filter(
            lambda x: len(x) >= 10
        )

        # compute average annual citation count
        avg_citation_count = compute_avg_citation_count(author_data)

        # compute sample publication count
        sample_publication_count = compute_sample_publication_count(author_data)

        # compute apf
        apf_data = compute_apf(author_data)

        # merge apf and avg_citation_count
        author_data = apf_data.merge(
            avg_citation_count, on=["author", "institution", "year"], how="left"
        )

        # merge with sample publication count
        author_data = author_data.merge(
            sample_publication_count, on=["author", "institution", "year"], how="left"
        )

        logger.info("Finished processing data from %s", key)

        yield {key: author_data}


def combine_lab_results(
    dict_loader: AbstractDataset,
) -> pd.DataFrame:
    """
    This node combines the lab determinants for each lab into a single DataFrame.

    Parameters:
        data (Dict[str, pd.DataFrame]): A dictionary containing the lab determinants
            for each lab.

    Returns:
        pd.DataFrame: A DataFrame containing the lab determinants for each lab.
    """
    data = []
    for i, (key, loader) in enumerate(dict_loader.items()):
        logger.info(
            "Loading data from %s. Number %s out of %s", key, i + 1, len(dict_loader)
        )
        data.append(loader())
    data = pd.concat(data, ignore_index=True)

    # change column "id" to "publication_count"
    data = data.rename(columns={"id": "publication_count"})

    # for duplicated author, institution, year triples, get the mean
    data = data.groupby(["author", "institution", "year"]).mean().reset_index()

    return data


def assign_lab_label(
    candidate_data: pd.DataFrame, quantile_val: float = 0.9
) -> pd.DataFrame:
    """
    Assigns lab labels to candidate data based on matching with ground truth data.

    Args:
        candidate_data (pd.DataFrame): The candidate data containing author information.
        quantile_val (float): The quantile value to use for filtering the candidate data.

    Returns:
        pd.DataFrame: The candidate data with lab labels assigned.

    """

    logger.info("Sorting candidate data")
    candidate_data = candidate_data.sort_values(by=["author", "year"]).reset_index(
        drop=True
    )

    logger.info("Filtering candidate data based on quantile 0.75 for publication count or cited by count at least one year before 2021")
    candidate_data["year"] = candidate_data["year"].astype(int)
    
    # filter data to years before 2021
    pre_2021_data = candidate_data[candidate_data["year"] < 2021]
    
    # calculate 0.75 quantiles for pre-2021 data
    pub_count_quantile_75 = pre_2021_data["publication_count"].quantile(0.75)
    cited_by_quantile_75 = pre_2021_data["cited_by_count"].quantile(0.75)
    
    # filter candidates who have at least one year before 2021 with either metric above 0.75 qt
    # keep all years (including post-2021) for qualifying candidates
    candidate_data = candidate_data.groupby(["author", "institution"]).filter(
        lambda group: (
            (group[group["year"] < 2021]["publication_count"] > pub_count_quantile_75).any() or
            (group[group["year"] < 2021]["cited_by_count"] > cited_by_quantile_75).any()
        )
    )

    logger.info(
        "Filtering candidate data based on appearing in three consecutive years"
    )
    candidate_data["year"] = candidate_data["year"].astype(int)
    candidate_data = candidate_data.groupby(["author", "institution"]).filter(
        appears_in_three_consecutive_years
    )

    logger.info("Calculating mean APF")
    candidate_data["mean_apf"] = candidate_data.groupby(["author", "institution"])[
        "apf"
    ].transform("mean")

    logger.info("Calculating 3-year average APF")
    candidate_data["apf_3yr_avg"] = (
        candidate_data.groupby(["author", "institution"])["apf"]
        .rolling(3)
        .mean()
        .reset_index(level=[0, 1], drop=True)
    )

    logger.info("Filling 1 and 2 year nans")
    candidate_data["apf_1yr_avg"] = (
        candidate_data.groupby(["author", "institution"])["apf"]
        .rolling(1)
        .mean()
        .reset_index(level=[0, 1], drop=True)
    )
    candidate_data["apf_2yr_avg"] = (
        candidate_data.groupby(["author", "institution"])["apf"]
        .rolling(2)
        .mean()
        .reset_index(level=[0, 1], drop=True)
    )

    # Use 1-year and 2-year averages to fill NaN in 'apf_3yr_avg'
    candidate_data["apf_3yr_avg"] = (
        candidate_data["apf_3yr_avg"]
        .fillna(candidate_data["apf_2yr_avg"])
        .fillna(candidate_data["apf_1yr_avg"])
    )

    # Drop temporary columns
    candidate_data = candidate_data.drop(columns=["apf_1yr_avg", "apf_2yr_avg"])

    logger.info("Normalising data")
    scaler = MinMaxScaler()

    def scale_group(group):
        group[["cited_by_count", "publication_count"]] = scaler.fit_transform(
            group[["cited_by_count", "publication_count"]]
        )
        return group

    candidate_data = (
        candidate_data.groupby("year").apply(scale_group).reset_index(drop=True)
    )

    # drop 2025 year
    candidate_data = candidate_data[candidate_data["year"] != 2025]

    # assign weights
    weights = {
        "apf": 0.3,
        "mean_apf": 0.1,
        "apf_3yr_avg": 0.2,
        "publication_count": 0.3,
        "cited_by_count": 0.1,
    }

    logger.info("Calculating score")
    candidate_data["score"] = (
        candidate_data["apf"] * weights["apf"]
        + candidate_data["mean_apf"] * weights["mean_apf"]
        + candidate_data["cited_by_count"] * weights["cited_by_count"]
        + candidate_data["publication_count"] * weights["publication_count"]
    )

    logger.info("Make final selection")
    quantile_res = candidate_data["score"].quantile(quantile_val)
    likely_pis = candidate_data.groupby(["author", "institution"]).filter(
        is_likely_pi, quantile=quantile_res
    )

    # Check if we have more than 15,000 unique authors and limit if necessary
    unique_authors = likely_pis["author"].nunique()
    if unique_authors > 15000:
        logger.info("Found %d unique authors, limiting to top 15,000 by author score", unique_authors)
        
        # calculate max score per author
        author_scores = likely_pis.groupby("author")["score"].max().reset_index()
        author_scores = author_scores.sort_values("score", ascending=False)
        
        # Select top 15,000 authors
        top_authors = author_scores.head(15000)["author"].tolist()
        
        # Filter likely_pis to only include these authors
        likely_pis = likely_pis[likely_pis["author"].isin(top_authors)]
        
        logger.info("Limited to %d unique authors", likely_pis["author"].nunique())

    return likely_pis


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

    filter_batches = [[from_publication_date, author_id] for author_id in author_ids]

    # slice to create parallel jobs that produce slices
    slices = [filter_batches[i : i + 150] for i in range(0, len(filter_batches), 150)]

    logger.info("Fetching papers for %d author batches", len(slices))

    for i, slice_ in enumerate(slices):

        logger.info("Processing batch number %d / %d", i + 1, len(slices))

        slice_results = Parallel(n_jobs=8, backend="loky", verbose=10)(
            delayed(collect_papers)(
                oa_ids=batch,
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
                        "publication_date",
                        "mesh_terms",
                        "cited_by_count",
                        "counts_by_year",
                        "authorships",
                        "concepts",
                        "topics",
                        "fwci",
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

                # break authorship nested dictionary jsons, create pipe-delimited string of authorship triplets
                df["authorships"] = df["authorships"].apply(
                    lambda x: (
                        "|".join([
                            f"{author['author']['id'].replace('https://openalex.org/', '')},{inst['id'].replace('https://openalex.org/', '') if inst else ''},{author['author_position']}"
                            for author in x
                            for inst in author["institutions"] or [{}]
                        ]) if x else None
                    )
                )

                df["primary_field"] = df["topics"].apply(
                    lambda x: (
                        x[0][5]
                        if isinstance(x, np.ndarray)
                        and len(x) > 0
                        and isinstance(x[0], np.ndarray)
                        and len(x[0]) > 0
                        else None
                    )
                )

                # create a list of topics
                try:
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
                except KeyError:
                    logger.warning("KeyError in topics for %s", df["id"].values[0])
                    df["topics"] = None

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

        # force new vars to float
        for col in ["fwci"]:
            # set any values with "" as NaN
            df[col] = df[col].replace("", np.nan).astype(float)

        yield {f"s{i}": df}


def create_candidates_map(
    alphafold_authors: List,
    ct_ai_authors: List,
    ct_pp_authors: List,
    ct_sb_authors: List,
    other_authors: List,
) -> pd.DataFrame:
    """
    Create a mapping between candidates and chain seeds.

    Args:
        alphafold_authors (List): The list of AlphaFold authors.
        ct_ai_authors (List): The list of CT AI authors.
        ct_pp_authors (List): The list of CT PP authors.
        ct_sb_authors (List): The list of CT SB authors.
        other_authors (List): The list of other authors.

    Returns:
        pd.DataFrame: A DataFrame containing the mapping between candidates and chain seeds.
    """
    # create a dataframe for each type of author
    alphafold_df = pd.DataFrame(alphafold_authors, columns=["author", "institution"])
    alphafold_df["seed"] = "alphafold"

    ct_ai_df = pd.DataFrame(ct_ai_authors, columns=["author", "institution"])
    ct_ai_df["seed"] = "ct_ai"

    ct_pp_df = pd.DataFrame(ct_pp_authors, columns=["author", "institution"])
    ct_pp_df["seed"] = "ct_pp"

    ct_sb_df = pd.DataFrame(ct_sb_authors, columns=["author", "institution"])
    ct_sb_df["seed"] = "ct_sb"

    other_df = pd.DataFrame(other_authors, columns=["author", "institution"])
    other_df["seed"] = "other"

    # concatenate all dataframes
    ct_df = pd.concat([ct_ai_df, ct_pp_df, ct_sb_df])

    # concatenate all dataframes
    candidate_map = pd.concat([alphafold_df, ct_df, other_df])

    # drop duplicates
    candidate_map = candidate_map.drop_duplicates(subset=["author", "institution"])

    # if seed contains ct_ai, go to ct_ai
    candidate_map["seed"] = candidate_map["seed"].apply(
        lambda x: "ct_ai" if "ct_ai" in x else x
    )

    # if seed contains ct_pp, go to ct_pp
    candidate_map["seed"] = candidate_map["seed"].apply(
        lambda x: "ct_pp" if "ct_pp" in x else x
    )

    # if seed contains ct_sb, go to ct_sb
    candidate_map["seed"] = candidate_map["seed"].apply(
        lambda x: "ct_sb" if "ct_sb" in x else x
    )

    return candidate_map


def get_institution_info(
    author_ids: pd.DataFrame,
) -> pd.DataFrame:
    """
    Fetches institution information for authors.

    Args:
        author_ids (pd.DataFrame): DataFrame containing author IDs and institution information.

    Returns:
        pd.DataFrame: DataFrame containing institution data.

    """
    institutions = author_ids["institution"].unique().tolist()

    logger.info("Fetching publications for %d authors", len(institutions))

    # slice to create parallel jobs that produce slices
    slices = [institutions[i : i + 150] for i in range(0, len(institutions), 150)]

    logger.info("Slicing authors into %d slices", len(slices))
    institution_data_list = Parallel(n_jobs=8, verbose=10)(
        delayed(fetch_institution_data)(institution)
        for slice_ in slices
        for institution in slice_
    )

    # filter out None values
    institution_data_list = [data for data in institution_data_list if data is not None]

    # convert the list of dictionaries to a DataFrame
    data = pd.DataFrame(institution_data_list)

    # merge back on author_ids
    data = author_ids[["author", "institution"]].merge(
        data, on="institution", how="left"
    )

    return data
