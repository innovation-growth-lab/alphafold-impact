"""
This is a boilerplate pipeline 'data_collection_ecr'
generated using Kedro 0.19.1
"""

import logging
from itertools import chain
from typing import Dict, Generator
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
from ..a_data_collection_oa.nodes import collect_papers  # pylint: disable=E0402


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
        break
    return authors_final


def fetch_ecr_outputs(
    authors: pd.DataFrame,
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
    ecr_authors = list(set(authors[authors["status"] == "ecr"]["candidate"].tolist()))

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

        slice_papers["authorships"] = slice_papers["authorships"].apply(
            lambda x: (
                [
                    author["author"].get("id", "").replace("https://openalex.org/", "")
                    for author in x
                ]
                if x
                else None
            )
        )

        # explode the authorships
        slice_papers = slice_papers.explode("authorships")

        # filter for authorships in the author list
        slice_papers = slice_papers[slice_papers["authorships"].isin(authors)]

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

        yield {f"s{i}": slice_papers}


def _result_transformations(data: pd.DataFrame) -> pd.DataFrame:
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

    # Extract the content of cited_by_percentile_year
    data[
        [
            "cited_by_percentile_year_min",
            "cited_by_percentile_year_max",
        ]
    ] = data.apply(
        lambda x: pd.Series(x["cited_by_percentile_year"]),
        axis=1,
        result_type="expand",
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
