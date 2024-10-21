"""
This is a boilerplate pipeline 'data_collection_ecr'
generated using Kedro 0.19.1
"""

import logging
from itertools import chain
from typing import Dict, Generator
from joblib import Parallel, delayed
from kedro.io import AbstractDataset
import pandas as pd
import numpy as np
from ..a_data_collection_oa.nodes import collect_papers  # pylint: disable=E0402
from ..e_data_output_publications.nodes import (  # pylint: disable=E0402
    get_patent_citations,
)


logger = logging.getLogger(__name__)


def _get_candidate_authors(data: pd.DataFrame) -> pd.DataFrame:
    """
    Get the candidate authors from the given DataFrame.

    Args:
        data (pd.DataFrame): The input DataFrame containing the data.

    Returns:
        pd.DataFrame: A DataFrame containing the candidate authors and their counts.
    """

    # get the authors
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
        .drop_duplicates(subset=["id", "author"])
        .reset_index(drop=True)
    )

    # drop A9999999999
    author_data = author_data[~(author_data["author"] == "A9999999999")]

    # create "foundational" or "applied" column
    author_data["depth"] = author_data["level"].apply(
        lambda x: "foundational" if x == "0" else "applied"
    )

    author_data = (
        author_data.groupby(["author", "institution", "depth", "source"])
        .size()
        .reset_index(name="counts")
    )

    return author_data


def get_unique_authors(
    publications_data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Extracts unique authors from publications data.

    Args:
        publications_data (pd.DataFrame): The input DataFrame containing the publications data.

    Returns:
        pd.DataFrame: A DataFrame with unique authors, their institutions, and their
                        aggregated counts from different sources.
    """

    logger.info("Getting unique authors")
    authors = _get_candidate_authors(
        publications_data[["id", "authorships", "source", "level"]]
    )

    # Find the institution with the highest count for each author
    top_institutions = authors.loc[
        authors.groupby(["author", "depth", "source"])["counts"].idxmax()
    ][["author", "depth", "source", "institution"]]
    authors = authors.drop("institution", axis=1).merge(
        top_institutions, on=["author", "depth", "source"]
    )

    # groupby author, institution and sum
    authors = (
        authors.groupby(["author", "institution", "depth", "source"])
        .sum()
        .reset_index()
    )

    # pivot the data
    authors = authors.pivot_table(
        index=["author", "institution", "depth"],
        columns="source",
        values="counts",
        fill_value=0,
    ).reset_index()

    # make source cols int
    for col in ["af", "ct", "other"]:
        authors[col] = authors[col].astype(int)

    return authors


def fetch_candidate_ecr_status(
    authors: pd.DataFrame,
    from_publication_date: str,
    to_publication_date: str,
    api_config: Dict[str, str],
) -> pd.DataFrame:
    """
    Fetches the status of candidate Early Career Researchers (ECR) based on their
    publication records. This function processes a list of authors, divides them
    into batches, and fetches their publication records within a specified date
    range using an API. It identifies authors who are not present in the fetched
    records, assuming them to be ECRs.

    Args:
        authors (pd.DataFrame): DataFrame containing author information with at
            least an 'author' column.
        from_publication_date (str): The start date for fetching publication records
            (format: 'YYYY-MM-DD').
        to_publication_date (str): The end date for fetching publication records
            (format: 'YYYY-MM-DD').
        api_config (Dict[str, str]): Configuration dictionary for the API, including
            parameters like 'perpage'.

    Returns:
        pd.DataFrame: A DataFrame containing the authors and their ECR status.
    """
    # create list of author ids
    author_ids = list(set(authors["author"].tolist()))

    # create batches of 50 authors
    author_batches = [
        "|".join(author_ids[i : i + 50]) for i in range(0, len(author_ids), 50)
    ]

    filter_batches = [
        [from_publication_date, to_publication_date, batch] for batch in author_batches
    ]

    # slice to create parallel jobs that produce slices
    slices = [filter_batches[i : i + 250] for i in range(0, len(filter_batches), 250)]

    # create authors final data
    authors_final = pd.DataFrame()

    logger.info("Fetching papers for %d author batches", len(slices))
    for i, slice_ in enumerate(slices):

        # create a unique list of the authors in the slice
        slice_candidates = list(chain.from_iterable([x[2].split("|") for x in slice_]))

        logger.info("Processing batch number %d / %d", i + 1, len(slices))

        slice_results = Parallel(n_jobs=8, backend="loky", verbose=10)(
            delayed(collect_papers)(
                oa_ids=batch,
                perpage=api_config["perpage"],
                filter_criteria=[
                    "from_publication_date",
                    "to_publication_date",
                    "author.id",
                ],
                slice_keys=True,
                eager_loading=True,
                skip_preprocess=True,
            )
            for batch in slice_
        )

        slice_papers = list(
            chain.from_iterable(
                paper["authorships"]
                for sublist in [v for d in slice_results for k, v in d.items()]
                for paper in sublist
            )
        )

        slice_papers = pd.DataFrame.from_records(slice_papers)

        slice_papers["author"] = slice_papers["author"].apply(
            lambda x: x.get("id", "").replace("https://openalex.org/", "")
        )

        slice_papers.drop_duplicates(subset=["author"], inplace=True)

        # Create a DataFrame from slice_candidates
        slice_candidates_df = pd.DataFrame(slice_candidates, columns=["candidate"])

        # Check if each candidate is in slice_papers and map the values
        slice_candidates_df["status"] = np.where(
            slice_candidates_df["candidate"].isin(slice_papers["author"]),
            "not_ecr",
            "ecr",
        )

        authors_final = pd.concat(
            [authors_final, slice_candidates_df], ignore_index=True
        )

    return authors_final


def fetch_author_outputs(
    authors: pd.DataFrame,
    ecr: bool,
    from_publication_date: str,
    api_config: Dict[str, str],
) -> Generator[Dict[str, pd.DataFrame], None, None]:
    """
    Fetches and processes author outputs from a specified publication date.

    Args:
        authors (pd.DataFrame): Author IDs to fetch outputs for.
        from_publication_date (str): The starting publication date to filter outputs.
        api_config (Dict[str, str]): Configuration dictionary for the API, including:
            - "perpage": Number of results per page.

    Yields:
        Dict[str, pd.DataFrame]: A dictionary where keys are slice identifiers and values
            are DataFrames containing the processed papers for each slice.
    """
    # create list of author ids
    if ecr:
        ecr_authors = list(
            set(authors[authors["status"] == "ecr"]["candidate"].tolist())
        )
    else:
        ecr_authors = list(
            set(authors[authors["status"] == "not_ecr"]["candidate"].tolist())
        )

    # create batches of 50 authors
    author_batches = [
        "|".join(ecr_authors[i : i + 50]) for i in range(0, len(ecr_authors), 50)
    ]

    filter_batches = [[from_publication_date, batch] for batch in author_batches]

    # slice to create parallel jobs that produce slices
    slices = [filter_batches[i : i + 250] for i in range(0, len(filter_batches), 250)]

    logger.info("Fetching papers for %d author batches", len(slices))
    for i, slice_ in enumerate(slices):

        # create a unique list of the authors in the slice
        authors = list(chain.from_iterable([x[1].split("|") for x in slice_]))

        logger.info("Processing batch number %d / %d", i + 1, len(slices))

        slice_results = Parallel(n_jobs=8, backend="loky", verbose=10)(
            delayed(collect_papers)(
                oa_ids=batch,
                perpage=api_config["perpage"],
                filter_criteria=[
                    "from_publication_date",
                    "author.id",
                ],
                slice_keys=True,
                eager_loading=True,
                skip_preprocess=True,
            )
            for batch in slice_
        )

        slice_papers = [
            paper
            for sublist in [v for d in slice_results for k, v in d.items()]
            for paper in sublist
        ]

        # match papers to the authors
        slice_papers = pd.DataFrame(slice_papers)

        # do necessary transformations
        slice_papers = _result_transformations(slice_papers)

        # explode the authorships
        slice_papers = slice_papers.explode("authorships")

        # separate tuples into columns
        slice_papers[["author", "institution", "author_position"]] = pd.DataFrame(
            slice_papers["authorships"].tolist(), index=slice_papers.index
        )

        # drop the authorships column
        slice_papers.drop(columns=["authorships"], inplace=True)

        # filter for authorships in the author list
        slice_papers = slice_papers[slice_papers["author"].isin(authors)]

        # normalise_citation_counts
        slice_papers = _normalise_citation_counts(slice_papers)

        # drop used variables
        slice_papers.drop(
            columns=[
                "counts_by_year",
                "cited_by_percentile_year",
                "citation_normalized_percentile",
            ],
            inplace=True,
        )

        # force fwci type to float
        slice_papers.replace({"fwci": ""}, np.nan, inplace=True)
        slice_papers["fwci"] = slice_papers["fwci"].astype(float)

        yield {f"s{i}": slice_papers}


def merge_author_data(
    data_loaders: Dict[str, AbstractDataset],
    candidate_authors: pd.DataFrame,
    authors: pd.DataFrame,
    ecr: bool,
    institutions: pd.DataFrame,
    patents_data: pd.DataFrame,
    pdb_submissions: pd.DataFrame,
    icite_data: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Merges ECR (Early Career Researcher) data with institution information, and
    other various data sources.

    Args:
        data_loaders (Dict[str, AbstractDataset]): A dictionary where keys are strings
            and values are data loader functions that return DataFrames.
        institutions (pd.DataFrame): A DataFrame containing institution information,
            which must include an 'institution' column.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames:
            - The first DataFrame contains the merged ECR data.
            - The second DataFrame contains the reduced ECR data for regression analysis.
    """

    institutions = _institution_preprocessing(institutions)
    candidate_authors = _candidates_preprocessing(candidate_authors)

    output_data = pd.DataFrame()
    for i, loader in enumerate(data_loaders.values()):
        logger.info(
            "Processing data loader %d / %d",
            i + 1,
            len(data_loaders),
        )
        data = loader()

        try:
            data.drop(
                columns=[
                    "concepts",
                    "mesh_terms",
                    "grants",
                    "topics",
                    "ids",
                    "cited_by_percentile_year_min",
                    "cited_by_percentile_year_max",
                    "cit_6",
                    "cit_7",
                ]
            )
        except:  # pylint: disable=bare-except
            pass

        data = data.merge(institutions, on="institution", how="left")
        data = data.merge(candidate_authors, on="author", how="left")

        # get patent citations
        data["doi"] = data["doi"].apply(
            lambda x: x.replace("https://doi.org/", "") if x is not None else None
        )
        data = get_patent_citations(data, patents_data)

        # get pdb activity metrics
        data = data.merge(
            pdb_submissions[["id", "resolution", "R_free"]], on="id", how="left"
        )

        # get icite_counts
        icite_outputs = _create_cc_counts(data[["id", "doi"]], icite_data)

        data = data.merge(icite_outputs, on="id", how="left")

        output_data = pd.concat([output_data, data], ignore_index=True)

    if ecr:
        # identify any authors that have snuck in from before 2020
        non_ecr_authors = output_data[
            output_data["publication_date"].apply(
                lambda x: pd.to_datetime(x) < pd.to_datetime("2020-01-01")
            )
        ]["author"].unique()

        # remove non-ecr authors
        output_data = output_data[~output_data["author"].isin(non_ecr_authors)]
    else:
        authors_list = authors[authors["status"] == "not_ecr"]["candidate"].tolist()
        non_ecr_authors = output_data[output_data["author"].isin(authors_list)][
            "author"
        ].unique()

        # remove non-ecr authors
        output_data = output_data[output_data["author"].isin(non_ecr_authors)]

    return output_data


def aggregate_to_quarterly(data: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregates the input data to the quarterly level.

    Args:
        data (pd.DataFrame): The input DataFrame containing the data
            to be aggregated.

    Returns:
        pd.DataFrame: The aggregated DataFrame with the data aggregated
            to the quarterly level.
    """
    data["publication_date"] = pd.to_datetime(data["publication_date"])
    data["quarter"] = data["publication_date"].dt.to_period("Q")
    for col in ["R_free", "resolution"]:
        data[col] = data[col].replace({"": np.nan}).astype("float")

    def safe_mode(series):
        mode = series.mode()
        return mode.iloc[0] if not mode.empty else np.nan

    return (
        data.groupby(["author", "quarter"])
        .agg(
            num_publications=("id", "size"),
            num_cited_by_count=("cited_by_count", "sum"),
            num_cit_0=("cit_0", "sum"),
            num_cit_1=("cit_1", "sum"),
            cited_by_count=("cited_by_count", "mean"),
            cit_0=("cit_0", "mean"),
            cit_1=("cit_1", "mean"),
            fwci=("fwci", "mean"),
            percentile_value=("citation_normalized_percentile_value", "mean"),
            patent_count=("patent_count", "sum"),
            patent_citation=("patent_citation", "sum"),
            ca_count=("ca_count", "sum"),
            resolution=("resolution", "mean"),
            R_free=("R_free", "mean"),
            num_publications_pdb=("R_free", lambda x: x.notna().sum()),
            institution=("institution", "first"),
            institution_cited_by_count=("institution_cited_by_count", "first"),
            country_code=("country_code", "first"),
            type=("type", "first"),
            depth=("depth", "first"),
            af=("af", "first"),
            ct=("ct", "first"),
            other=("other", "first"),
            primary_field=("primary_field", safe_mode),
            author_position=("author_position", safe_mode),
        )
        .reset_index(),
    )


def _institution_preprocessing(institutions: pd.DataFrame) -> pd.DataFrame:
    institutions = institutions.drop(columns=["author"]).drop_duplicates(
        subset=["institution"]
    )
    institutions.rename(
        columns={"cited_by_count": "institution_cited_by_count"}, inplace=True
    )

    return institutions


def _candidates_preprocessing(candidate_authors: pd.DataFrame) -> pd.DataFrame:
    # sort candidate authors by af, ct, other
    depth_summed_data = (
        candidate_authors.groupby(["author", "depth"])
        .agg({"af": "sum", "ct": "sum", "other": "sum"})
        .reset_index()
    )

    depth_summed_data["total"] = (
        depth_summed_data["af"] + depth_summed_data["ct"] + depth_summed_data["other"]
    )

    # Get the index of the row with the maximum total for each author
    idx = depth_summed_data.groupby("author")["total"].idxmax()

    max_total_data = depth_summed_data.loc[idx]

    summed_data = (
        candidate_authors.groupby("author")
        .agg({"af": "sum", "ct": "sum", "other": "sum"})
        .reset_index()
    )

    # merge the summed data with the max total data
    max_total_data = max_total_data[["author", "depth"]].merge(
        summed_data, on="author", how="left"
    )

    return max_total_data


def _get_mesh_data(data: pd.DataFrame, mesh_terms: pd.DataFrame) -> pd.DataFrame:

    mesh_terms["term_group"] = mesh_terms["tree_number"].apply(
        lambda x: str(x)[:1] if x is not None else None
    )
    mesh_tag_class_dict = mesh_terms.set_index("DUI")["term_group"].to_dict()

    data["major_mesh_terms"] = data["mesh_terms"].apply(
        lambda x: (
            _extract_mesh_terms(x, True, mesh_tag_class_dict)
            if isinstance(x, np.ndarray)
            else []
        )
    )
    data["minor_mesh_terms"] = data["mesh_terms"].apply(
        lambda x: (
            _extract_mesh_terms(x, False, mesh_tag_class_dict)
            if isinstance(x, np.ndarray)
            else []
        )
    )

    return data


def _extract_mesh_terms(mesh_terms, is_major, mesh_dict):
    terms = []
    for term in mesh_terms:
        if term["is_major_topic"] == is_major:
            descriptor_ui = term["descriptor_ui"]
            descriptor_name = term["descriptor_name"]
            mesh_class = mesh_dict.get(descriptor_ui, None)
            methods = term["qualifier_name"]
            terms.append([descriptor_name, mesh_class, methods])
    return terms


def _create_cc_counts(data, icite_data):
    """
    Creates a count column in the given data DataFrame based on the 'cited_by_clin' column
        in the icite_data DataFrame.

    Args:
        data (pandas.DataFrame): The input data DataFrame.
        icite_data (pandas.DataFrame): The icite_data DataFrame containing the
            'cited_by_clin' column.

    Returns:
        pandas.DataFrame: The updated data DataFrame with a new 'count' column.
    """
    # change pmid, doi to str
    icite_data["doi"] = icite_data["doi"].astype(str)

    # merge on 'doi'
    data_doi = data.merge(icite_data[["doi", "cited_by_clin"]], how="left", on="doi")
    data_doi = data_doi.drop_duplicates(subset=["id"])

    logger.info("Exploding cited_by_clin")
    data_doi = data_doi[data_doi["cited_by_clin"].astype(str) != "nan"]
    data_doi = data_doi[data_doi["cited_by_clin"].notnull()]
    data_doi["cited_by_clin"] = data_doi["cited_by_clin"].apply(lambda x: x.split(" "))

    # create count column
    data_doi["ca_count"] = data_doi["cited_by_clin"].apply(len)

    # merge back with data, fill with 0 for missing count
    data = data.merge(data_doi[["id", "ca_count"]], on="id", how="left")
    data["ca_count"] = data["ca_count"].fillna(0)

    return data[["id", "ca_count"]]


def _result_transformations(data: pd.DataFrame) -> pd.DataFrame:
    """
    Transforms the input DataFrame by extracting and restructuring specific fields.
    The function performs the following transformations:
        1. Extracts and restructures the 'authorships' field.
        2. Creates a list of topics from the 'topics' field.
        3. Extracts concepts from the 'concepts' field.
        4. Extracts the content of 'citation_normalized_percentile' into separate columns.
        5. Extracts the content of 'cited_by_percentile_year' into separate columns.

    Args:
        data (pd.DataFrame): The input DataFrame containing the data to be transformed.

    Returns:
        pd.DataFrame: The transformed DataFrame with the extracted and restructured fields.
    """
    # extract the content of authorships
    data["authorships"] = data["authorships"].apply(
        lambda x: (
            [
                (
                    (
                        author["author"]["id"].replace("https://openalex.org/", ""),
                        inst["id"].replace("https://openalex.org/", ""),
                        author["author_position"],
                    )
                    if author["institutions"]
                    else [
                        author["author"]["id"].replace("https://openalex.org/", ""),
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
    data["topics"] = data["topics"].apply(
        lambda x: (
            [
                (
                    topic.get("id", "").replace("https://openalex.org/", ""),
                    topic.get("display_name", ""),
                    topic.get("subfield", {})
                    .get("id", "")
                    .replace("https://openalex.org/", ""),
                    topic.get("subfield", {}).get("display_name", ""),
                    topic.get("field", {})
                    .get("id", "")
                    .replace("https://openalex.org/", ""),
                    topic.get("field", {}).get("display_name", ""),
                    topic.get("domain", {})
                    .get("id", "")
                    .replace("https://openalex.org/", ""),
                    topic.get("domain", {}).get("display_name", ""),
                )
                for topic in x
            ]
            if x
            else None
        )
    )

    # primary field
    data["primary_field"] = data["topics"].apply(
        lambda x: (
            x[0][5]
            if isinstance(x, np.ndarray)
            and len(x) > 0
            and isinstance(x[0], np.ndarray)
            and len(x[0]) > 0
            else None
        )
    )

    # extract concepts
    data["concepts"] = data["concepts"].apply(
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

    # Extract the content of citation_normalized_percentile
    data[
        [
            "citation_normalized_percentile_value",
            "citation_normalized_percentile_is_in_top_1_percent",
            "citation_normalized_percentile_is_in_top_10_percent",
        ]
    ] = data.apply(
        lambda x: (pd.Series(x["citation_normalized_percentile"])),
        axis=1,
        result_type="expand",
    )

    data.drop(
        columns=[
            "concepts",
            "mesh_terms",
            "grants",
            "topics",
            "ids",
            "cited_by_percentile_year_min",
            "cited_by_percentile_year_max",
        ]
    )

    return data


def _normalise_citation_counts(data: pd.DataFrame) -> pd.DataFrame:
    """
    Normalises the citation counts in the given DataFrame.

    Args:
        data (pd.DataFrame): The input DataFrame containing the data.

    Returns:
        pd.DataFrame: A DataFrame with normalised citation counts.
    """

    data = data.drop_duplicates(subset=["id"])

    # explode the 'counts_by_year' column into multiple rows
    data_exploded = data.explode("counts_by_year")

    data_exploded.reset_index(drop=True, inplace=True)

    # separate the dictionary into two columns: 'year' and 'cited_by_count'
    data_exploded[["year", "cited_by_count"]] = pd.json_normalize(
        data_exploded["counts_by_year"]
    )

    # convert the 'publication_date' column to datetime and extract the year
    data_exploded["publication_date"] = pd.to_datetime(
        data_exploded["publication_date"]
    ).dt.year

    # calculate the 't' value
    data_exploded["t"] = data_exploded["year"] - data_exploded["publication_date"]

    # drop the rows with NaN values in the 't' column
    data_exploded.dropna(subset=["t"], inplace=True)

    data_exploded["t"] = data_exploded["t"].astype(int)

    # pivot the DataFrame to get the 't' columns
    t_columns = data_exploded.pivot_table(
        index="id", columns="t", values="cited_by_count", aggfunc="sum"
    )

    # remove negative columns
    t_columns = t_columns.loc[:, (t_columns.columns >= 0)]

    # Convert columns to Int64 type to accept NAs
    t_columns = t_columns.astype(pd.Int64Dtype())

    # renaming
    t_columns.columns = ["cit_" + str(int(i)) for i in t_columns.columns]

    # fill remaining NaN values in 't0' with 0
    t_columns["cit_0"] = t_columns["cit_0"].fillna(0)

    # merge the 't' DataFrame back to the original DataFrame
    data = data.merge(t_columns, left_on="id", right_index=True, how="left")

    return data
