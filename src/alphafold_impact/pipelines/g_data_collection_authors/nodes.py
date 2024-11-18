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

pd.set_option("future.no_silent_downcasting", True)
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
        .drop(columns=["authorships"])
    )

    # drop A9999999999
    author_data = author_data[~(author_data["author"] == "A9999999999")]
    author_data["level"] = author_data["level"].astype(int)

    # create "foundational" or "applied" column
    author_data["depth"] = author_data["level"].apply(
        lambda x: "foundational" if x == 0 else "applied"
    )

    # create quarter column
    author_data["publication_date"] = pd.to_datetime(author_data["publication_date"])
    author_data["quarter"] = author_data["publication_date"].dt.to_period("Q")

    author_data_out = (
        author_data.groupby(["author", "institution", "depth", "source", "quarter"])
        .size()
        .reset_index(name="counts")
    )

    logger.info("Getting chain labels")
    sort_order = {
        "strong": 0,
        "partial_strong": 1,
        "mixed": 2,
        "weak": 3,
        "partial_weak": 4,
        "unknown": 5,
        "no_data": 6,
    }
    author_data["sort_order"] = author_data.apply(
        lambda row: -1 if row["level"] == -1 else sort_order.get(row["chain_label"], 7),
        axis=1,
    )

    author_labels = (
        author_data[
            ["sort_order", "quarter", "author", "chain_label", "level", "source"]
        ]
        .sort_values(["level", "sort_order"])
        .groupby(["author", "quarter", "source"])
        .first()
        .reset_index()
        .drop(columns=["sort_order", "level"])
    )

    # keep rows with strong or partial storng chain labels
    author_labels = author_labels[
        author_labels["chain_label"].isin(["strong", "partial_strong"])
    ]

    author_labels["strong"] = author_labels["source"]

    return author_data_out, author_labels


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
    authors, strength_labels = _get_candidate_authors(
        publications_data[
            ["id", "authorships", "source", "level", "publication_date", "chain_label"]
        ]
    )

    logger.info("Getting unique institutions")
    top_institutions = authors.loc[
        authors.groupby(["author", "depth", "source"])["counts"].idxmax()
    ][["author", "depth", "source", "institution"]]
    authors = authors.drop("institution", axis=1).merge(
        top_institutions, on=["author", "depth", "source"]
    )

    logger.info("Cumulative sum of counts")
    authors = (
        authors.groupby(["author", "depth", "source", "quarter"])["counts"]
        .sum()
        .reset_index()
    )
    authors = authors.sort_values(by=["author", "depth", "source", "quarter"])
    authors["cumulative_counts"] = authors.groupby(["author", "depth", "source"])[
        "counts"
    ].cumsum()

    logger.info("Pivoting table")
    authors = authors.pivot_table(
        index=["author", "depth", "quarter"],
        columns="source",
        values="cumulative_counts",
        fill_value=0,
    ).reset_index()

    # make source cols int
    for col in ["af", "ct_ai", "ct_noai", "other"]:
        authors[col] = authors[col].astype(int)

    logger.info("Pivoting labels")
    strength_labels_wide = _pivot_labels(strength_labels)

    authors = authors.merge(strength_labels_wide, on=["author", "quarter"], how="left")

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
    slices = [filter_batches[i : i + 100] for i in range(0, len(filter_batches), 100)]

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
    mesh_terms: pd.DataFrame,
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
    logger.info("Preparing labels for end-of-loop merge")
    label_cols = [col for col in candidate_authors.columns if col.startswith("strong")]
    author_labels = candidate_authors[
        ["author", "quarter"] + label_cols
    ].drop_duplicates(subset=["author", "quarter"])

    # sort by author and quarter
    author_labels = author_labels.sort_values(by=["author", "quarter"])

    # ffill strong labels
    author_labels[
        [col for col in author_labels.columns if col.startswith("strong")]
    ] = (author_labels[["author"] + label_cols].groupby("author").ffill().fillna(0))

    for col in label_cols:
        author_labels[col] = author_labels[col].astype(int)

    candidate_authors = candidate_authors.drop(
        columns=[col for col in candidate_authors.columns if col.startswith("strong")]
    )

    institutions = _institution_preprocessing(institutions)
    candidate_authors = _candidates_preprocessing(candidate_authors)

    mesh_terms_dict = mesh_terms.set_index("DUI")["term_group"].to_dict()

    output_data = pd.DataFrame()
    for i, loader in enumerate(data_loaders.values()):
        logger.info(
            "Processing data loader %d / %d",
            i + 1,
            len(data_loaders),
        )
        data = loader()

        logger.info("Getting high PDB authors")
        data = _define_high_pdb_authors(data, pdb_submissions)

        # get first author data
        data = data[data["author_position"] == "first"]

        data["publication_date"] = pd.to_datetime(data["publication_date"])
        data["quarter"] = data["publication_date"].dt.to_period("Q")

        # HACK: ecr data has no primary_field during collection
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

        logger.info("Calculating field share")
        data = calculate_field_share(data)

        logger.info("Calculating mesh balance")
        data = calculate_mesh_balance(data, mesh_terms_dict)

        logger.info("Collecting COVID references")
        data = collect_covid_references(data)

        data.drop(
            columns=[
                "concepts",
                "mesh_terms",
                "grants",
                "topics",
                "ids",
            ],
            inplace=True,
        )

        if ecr:
            # identify any authors that have snuck in from before 2020
            non_ecr_authors = data[
                data["publication_date"].apply(
                    lambda x: pd.to_datetime(x) < pd.to_datetime("2020-01-01")
                )
            ]["author"].unique()

            # remove non-ecr authors
            data = data[~data["author"].isin(non_ecr_authors)]
        else:
            authors_list = authors[authors["status"] == "not_ecr"]["candidate"].tolist()
            non_ecr_authors = data[data["author"].isin(authors_list)]["author"].unique()

            # remove ecr authors
            data = data[data["author"].isin(non_ecr_authors)]

        # merge author and institutions metadata
        data = data.merge(institutions, on="institution", how="left")
        data = data.merge(candidate_authors, on=["author", "quarter"], how="left")

        # sort by quarter, bfill and ffill for missing af. ct_ai, ct_noai, other
        data = data.sort_values(by=["author", "quarter"])
        data[["af", "ct_ai", "ct_noai", "other"]] = (
            data.groupby("author")[["af", "ct_ai", "ct_noai", "other"]]
            .ffill()
            .fillna(0)
            .astype(int)
        )

        # fill depth with whatever value groupby author is not NAN
        data["depth"] = data.groupby("author")["depth"].transform(
            lambda x: x.ffill().bfill()
        )

        # drop if depth still nan
        data = data.dropna(subset=["depth"])

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

        # concatenate the data
        output_data = pd.concat([output_data, data], ignore_index=True)

    # merge the unique authors with the output data
    output_data = output_data.merge(author_labels, on=["author", "quarter"], how="left")

    return output_data


def _institution_preprocessing(institutions: pd.DataFrame) -> pd.DataFrame:
    institutions = institutions.drop(columns=["author"]).drop_duplicates(
        subset=["institution"]
    )
    institutions.rename(
        columns={"cited_by_count": "institution_cited_by_count"}, inplace=True
    )

    return institutions


def _candidates_preprocessing(candidate_authors: pd.DataFrame) -> pd.DataFrame:
    """
    Create a DataFrame with the relative sum for each author in each depth.
    Assign the depth to the author based on the relative sum where the author
    is represented most often.
    """
    # Calculate the total sum for each depth
    total_depth_sums = (
        candidate_authors.groupby("depth")
        .agg({"af": "sum", "ct_ai": "sum", "ct_noai": "sum", "other": "sum"})
        .reset_index()
    )

    total_depth_sums["total"] = (
        total_depth_sums["af"]
        + total_depth_sums["ct_ai"]
        + total_depth_sums["ct_noai"]
        + total_depth_sums["other"]
    )

    # Calculate the relative sum for each author in each depth
    author_depth_sums = (
        candidate_authors.groupby(["author", "depth"])
        .agg({"af": "sum", "ct_ai": "sum", "ct_noai": "sum", "other": "sum"})
        .reset_index()
    )

    author_depth_sums["total"] = (
        author_depth_sums["af"]
        + author_depth_sums["ct_ai"]
        + author_depth_sums["ct_noai"]
        + author_depth_sums["other"]
    )

    # Merge the total depth sums with the author depth sums
    relative_sums = author_depth_sums.merge(
        total_depth_sums[["depth", "total"]], on="depth", suffixes=("", "_depth")
    )

    # Calculate the relative relevance
    relative_sums["relative_relevance"] = (
        relative_sums["total"] / relative_sums["total_depth"]
    )

    # Get the index of the row with the maximum relative relevance for each author
    idx = relative_sums.groupby("author")["relative_relevance"].idxmax()

    max_relative_data = relative_sums.loc[idx]

    # Group by author and quarter to sum the values
    summed_data = (
        candidate_authors.groupby(["author", "quarter"])
        .agg({"af": "sum", "ct_ai": "sum", "ct_noai": "sum", "other": "sum"})
        .reset_index()
    )

    # Merge the summed data with the max relative data
    max_relative_data = max_relative_data[["author", "depth"]].merge(
        summed_data, on="author", how="left"
    )

    return max_relative_data


def _define_high_pdb_authors(
    data: pd.DataFrame, pdb_submissions: pd.DataFrame
) -> pd.DataFrame:
    """
    Define high PDB authors based on the publication counts in the
    fourth quantile.

    Args:
        data (pd.DataFrame): The input DataFrame containing the data.
        pdb_submissions (pd.DataFrame): The DataFrame containing the PDB submissions.

    Returns:
        pd.DataFrame: The updated DataFrame with the 'high_pdb' column.
    """

    data_db = data[["id", "author"]].merge(pdb_submissions, on="id", how="inner")

    data_db["pdb_submission"] = True

    # Filter data_db to include only publications before 2021
    data_db_pre_2021 = data_db[data_db["publication_date"] < "2021-01-01"]

    # Group by author and count the number of publications for each author
    author_pub_counts = data_db_pre_2021.groupby("author")["id"].count().reset_index()
    author_pub_counts = author_pub_counts.rename(columns={"id": "pub_count"})

    # Calculate the 75th percentile (fourth quantile) of the publication counts
    quantile_75 = author_pub_counts["pub_count"].quantile(0.75)

    # Create a high_pdb variable for authors with publication counts in the fourth quantile
    author_pub_counts["high_pdb"] = author_pub_counts["pub_count"] >= quantile_75

    # Merge this information back into the data DataFrame
    data = data.merge(
        author_pub_counts[["author", "high_pdb"]], on="author", how="left"
    )

    # Merge data_db to include pdb_submission column
    data = data.merge(data_db[["id", "pdb_submission"]], on="id", how="left")
    data.fillna({"pdb_submission": False}, inplace=True)

    return data


def calculate_field_share(data):
    """
    Calculate the share of each subfield for each author in a DataFrame.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data.

    Returns:
        pandas.DataFrame: The resulting DataFrame with the share of each subfield for each author.
    """
    df = data.copy()
    df["fields"] = df["topics"].apply(
        lambda x: [y[5] for y in x] if x is not None and len(x) > 0 else []
    )
    # explode the subfields column into multiple rows
    df = df.explode("fields")
    # group by author and subfields and calculate the count
    df = df.groupby(["author", "quarter", "fields"]).size().reset_index(name="count")
    # calculate the total count for each author
    total_count = df.groupby("author")["count"].sum()
    # calculate the share of each subfield
    df["share"] = df.apply(
        lambda row: row["count"] / total_count[row["author"]], axis=1
    )
    # pivot the DataFrame to get one column for each subfield
    df = df.pivot(index=["author", "quarter"], columns="fields", values="share")

    # reset index
    df.reset_index(inplace=True)

    # fill NaN values with 0
    df = df.fillna(0)

    # only keep first 10 characters of the subfield
    df.columns = [
        col[:10] if col != "author" and col != "quarter" else col for col in df.columns
    ]

    # change column names to be camel case, and prefix with "field_"
    df.columns = [
        (
            "field_" + col.lower().replace(" ", "_")
            if col != "author" and col != "quarter"
            else col
        )
        for col in df.columns
    ]

    # try to drop "field_" column
    if "field_" in df.columns:
        df = df.drop(columns=["field_"])

    # merge with data on author and time
    data = data.merge(df, on=["author", "quarter"], how="left")

    data.drop(columns=["topics"], inplace=True)

    return data


def calculate_mesh_balance(
    data: pd.DataFrame, mesh_terms: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate the balance of mesh terms for each author over time.

        Args:
        data (pd.DataFrame): A DataFrame containing author data with columns 'author',
          'time', and 'mesh_terms'.
        mesh_terms (pd.DataFrame): A DataFrame containing mesh terms and their
          corresponding groups.
    Returns:
        pd.DataFrame: A DataFrame with the original data and additional columns for
          each mesh term's share, prefixed with 'mesh_'.
    """

    df = data.copy()

    # transform mesh_terms by mapping to term_group
    df["mesh_terms"] = df["mesh_terms"].apply(
        lambda x: [y[0] for y in x] if x is not None else []
    )
    df["mesh_terms"] = df["mesh_terms"].apply(
        lambda x: [mesh_terms.get(y, "") for y in x] if x is not None else []
    )

    # explode the mesh_terms column into multiple rows
    df = df.explode("mesh_terms")

    # group by author and mesh_terms and calculate the count
    df = (
        df.groupby(["author", "quarter", "mesh_terms"]).size().reset_index(name="count")
    )

    # calculate the total count for each author
    total_count = df.groupby("author")["count"].sum()

    # calculate the share of each mesh term
    df["share"] = df.apply(
        lambda row: row["count"] / total_count[row["author"]], axis=1
    )

    # pivot the DataFrame to get one column for each mesh term
    df = df.pivot(index=["author", "quarter"], columns="mesh_terms", values="share")

    # reset index
    df.reset_index(inplace=True)

    # change column names to be camel case, and prefix with "mesh_"
    df.columns = [
        "mesh_" + col if col != "author" and col != "quarter" else col
        for col in df.columns
    ]

    # fill NaN values with 0
    df = df.fillna(0)

    # try to drop "mesh_" column
    if "mesh_" in df.columns:
        df = df.drop(columns=["mesh_"])

    # merge with data on author and time
    data = data.merge(df, on=["author", "quarter"], how="left")

    data.drop(columns=["mesh_terms"], inplace=True)

    return data


def collect_covid_references(data: pd.DataFrame) -> pd.DataFrame:
    """
    Collects and calculates the share of COVID-19 related concepts for each principal
     investigator (author) in the provided DataFrame for the months of September to
     December 2020.

    Args:
      data (pd.DataFrame): A DataFrame containing at least the following columns:
    Returns:
      pd.DataFrame: The original DataFrame with an additional column "covid_share_2020"
    """

    data_2020 = data[data["quarter"].isin(["2020Q1", "2020Q2", "2020Q3", "2020Q4"])]

    data_2020["covid_2020"] = data_2020["concepts"].apply(
        lambda x: (any(element == "C524204448" for element in x))
    )

    # Group the filtered dataframe by author and calculate the share of COVID concepts
    covid_share_2020 = (
        data_2020.groupby("author")["covid_2020"].sum()
        / data_2020.groupby("author")["covid_2020"].count()
    )

    # Convert the Series to a DataFrame
    covid_share_2020 = covid_share_2020.reset_index()

    # Rename the column
    covid_share_2020.columns = ["author", "covid_share_2020"]

    # Create a dictionary from the DataFrame
    covid_share_dict = dict(
        zip(covid_share_2020["author"], covid_share_2020["covid_share_2020"])
    )

    data["covid_share_2020"] = data["author"].map(covid_share_dict)

    data.drop(columns=["concepts"], inplace=True)

    return data


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
    for col in ["R_free", "resolution"]:
        data[col] = data[col].replace({"": np.nan}).astype("float")

    def safe_mode(series):
        mode = series.mode()
        return mode.iloc[0] if not mode.empty else np.nan

    # aggregation dictionaryy
    agg_dict = {
        "num_publications": ("id", "size"),
        "num_cited_by_count": ("cited_by_count", "sum"),
        "num_cit_0": ("cit_0", "sum"),
        "num_cit_1": ("cit_1", "sum"),
        "cited_by_count": ("cited_by_count", "mean"),
        "cit_0": ("cit_0", "mean"),
        "cit_1": ("cit_1", "mean"),
        "fwci": ("fwci", "mean"),
        "percentile_value": ("citation_normalized_percentile_value", "mean"),
        "patent_count": ("patent_count", "sum"),
        "patent_citation": ("patent_citation", "sum"),
        "ca_count": ("ca_count", "sum"),
        "resolution": ("resolution", "mean"),
        "R_free": ("R_free", "mean"),
        "num_publications_pdb": ("pdb_submission", "sum"),
        "institution": ("institution", "first"),
        "institution_cited_by_count": ("institution_cited_by_count", "first"),
        "country_code": ("country_code", "first"),
        "type": ("type", "first"),
        "depth": ("depth", "first"),
        "af": ("af", "first"),
        "ct_ai": ("ct_ai", "first"),
        "ct_noai": ("ct_noai", "first"),
        "other": ("other", "first"),
        "strong_af": ("strong_af", safe_mode),
        "strong_ct_ai": ("strong_ct_ai", safe_mode),
        "strong_ct_noai": ("strong_ct_noai", safe_mode),
        "strong_other": ("strong_other", safe_mode),
        "primary_field": ("primary_field", safe_mode),
        "author_position": ("author_position", safe_mode),
        "high_pdb": ("high_pdb", "first"),
    }

    # add fields, mesh, covid
    additional_columns = [
        col
        for col in data.columns
        if col.startswith("field_") or col.startswith("mesh_")
    ]
    additional_columns.append("covid_share_2020")

    for col in additional_columns:
        agg_dict[col] = (col, "first")

    # group by author and quarter, aggregate calcs
    output_data = data.groupby(["author", "quarter"]).agg(**agg_dict).reset_index()

    return output_data


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
            # "concepts",
            # "mesh_terms",
            "grants",
            # "topics",
            "ids",
        ],
        inplace=True,
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


def _pivot_labels(data):
    strength_labels_wide = data.pivot_table(
        index=["author", "quarter"], columns="source", values="strong", aggfunc="first"
    )

    # relabel
    strength_labels_wide.columns = [
        f"strong_{col}" for col in strength_labels_wide.columns
    ]
    cols = [col for col in strength_labels_wide.columns]

    # sort
    strength_labels_wide = strength_labels_wide.reset_index().sort_values(
        by=["author", "quarter"]
    )

    for col in cols:
        strength_labels_wide[col] = np.where(strength_labels_wide[col].notnull(), 1, 0)

    # get the cumulative sum of counts
    strength_labels_wide[cols] = strength_labels_wide.groupby("author")[cols].cumsum()

    # forward fill
    strength_labels_wide[cols] = strength_labels_wide.groupby("author")[cols].ffill()

    return strength_labels_wide
