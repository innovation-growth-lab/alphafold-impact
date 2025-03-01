"""
This module contains functions for data collection and processing of author information
from various sources. The functions include extracting unique authors, fetching candidate
ECR status, fetching author outputs, merging author data, and various helper functions
for data transformations and calculations.

Functions:
    _get_candidate_authors(data: pd.DataFrame) -> pd.DataFrame:

    get_unique_authors(publications_data: pd.DataFrame) -> pd.DataFrame:

    fetch_candidate_ecr_status(authors: pd.DataFrame, from_publication_date: str,
        to_publication_date: str, api_config: Dict[str, str]) -> pd.DataFrame:
        publication records.

    fetch_author_outputs(authors: pd.DataFrame, ecr: bool, from_publication_date: str,
        api_config: Dict[str, str]) -> Generator[Dict[str, pd.DataFrame], None, None]:

    merge_author_data(data_loaders: Dict[str, AbstractDataset], candidate_authors: pd.DataFrame,
        authors: pd.DataFrame, ecr: bool, institutions: pd.DataFrame,
        patents_data: pd.DataFrame, pdb_submissions: pd.DataFrame,
        icite_data: pd.DataFrame, mesh_terms: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:

    _institution_preprocessing(institutions: pd.DataFrame) -> pd.DataFrame:
        Preprocesses the institution data.

    _candidates_preprocessing(candidate_authors: pd.DataFrame) -> pd.DataFrame:
        Preprocesses the candidate authors data.

    define_high_pdb_authors(data: pd.DataFrame, pdb_submissions: pd.DataFrame) -> pd.DataFrame:
        Define high PDB authors based on the publication counts in the fourth quantile.

    calculate_field_share(data: pd.DataFrame) -> pd.DataFrame:

    calculate_mesh_balance(data: pd.DataFrame, mesh_terms: pd.DataFrame) -> pd.DataFrame:

    collect_covid_references(data: pd.DataFrame) -> pd.DataFrame:

    aggregate_to_quarterly(data: pd.DataFrame) -> pd.DataFrame:

    create_cc_counts(data: pd.DataFrame, icite_data: pd.DataFrame) -> pd.DataFrame:

    _result_transformations(data: pd.DataFrame) -> pd.DataFrame:

    _normalise_citation_counts(data: pd.DataFrame) -> pd.DataFrame:

    _pivot_labels(data: pd.DataFrame) -> pd.DataFrame:
        Pivots the strength labels to a wide format and performs cumulative sum and forward fill.

"""

import logging
from itertools import chain
from typing import Dict, Generator
from joblib import Parallel, delayed
from kedro.io import AbstractDataset
import pandas as pd
import numpy as np
from ..a_data_collection_oa.nodes import collect_papers  # pylint: disable=E0402

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

    # create quarter column
    author_data["publication_date"] = pd.to_datetime(author_data["publication_date"])
    author_data["quarter"] = author_data["publication_date"].dt.to_period("Q")

    author_data_out = (
        author_data.groupby(["author", "institution", "source", "quarter"])
        .size()
        .reset_index(name="counts")
    )

    # create grouped chain_labels
    labels = {
        "strong": ["strong", "partial_strong"],
        "weak": ["weak", "partial_weak"],
        "mixed": ["mixed"],
        "unknown": ["unknown", "no_data"],
    }

    # create author labels
    author_data["grouped_chain_label"] = author_data["chain_label"].apply(
        lambda x: next((k for k, v in labels.items() if x in v), "unknown")
    )

    author_labels = (
        author_data.groupby(["author", "quarter", "grouped_chain_label", "source"])
        .size()
        .reset_index(name="counts")
    )

    # pivot the author labels
    author_labels = author_labels.pivot_table(
        index=["author", "quarter"],
        columns=["source", "grouped_chain_label"],
        values="counts",
        fill_value=0,
    ).reset_index()

    # join the two-level columns into one
    author_labels.columns = [
        f"{col[0]}_{col[1]}" if isinstance(col, tuple) else col
        for col in author_labels.columns
    ]

    # rename author_ and quarter_ to author and quarter
    author_labels.rename(
        columns={"author_": "author", "quarter_": "quarter"}, inplace=True
    )

    # make all columns but author_ and quarter_ int
    author_labels = author_labels.astype(
        {
            col: "int"
            for col in author_labels.columns
            if col not in ["author", "quarter"]
        }
    )

    # create four columns that sum the other counts on row and substract the count of _unknown_
    author_labels["af_with_intent"] = (
        author_labels["af_strong"]
        + author_labels["af_weak"]
        + author_labels["af_mixed"]
    )
    author_labels["ct_ai_with_intent"] = (
        author_labels["ct_ai_strong"]
        + author_labels["ct_ai_weak"]
        + author_labels["ct_ai_mixed"]
    )
    author_labels["ct_noai_with_intent"] = (
        author_labels["ct_noai_strong"]
        + author_labels["ct_noai_weak"]
        + author_labels["ct_noai_mixed"]
    )
    author_labels["other_with_intent"] = (
        author_labels["other_strong"]
        + author_labels["other_weak"]
        + author_labels["other_mixed"]
    )

    count_cols = [
        col for col in author_labels.columns if col not in ["author", "quarter"]
    ]
    author_labels[count_cols] = author_labels.groupby("author")[count_cols].cumsum()

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
        authors.groupby(["author", "source"])["counts"].idxmax()
    ][["author", "source", "institution"]]
    authors = authors.drop("institution", axis=1).merge(
        top_institutions, on=["author", "source"]
    )

    logger.info("Cumulative sum of counts")
    authors = (
        authors.groupby(["author", "source", "quarter"])["counts"].sum().reset_index()
    )
    authors = authors.sort_values(by=["author", "source", "quarter"])
    authors["cumulative_counts"] = authors.groupby(["author", "source"])[
        "counts"
    ].cumsum()

    logger.info("Pivoting table")
    authors = authors.pivot_table(
        index=["author", "quarter"],
        columns="source",
        values="cumulative_counts",
        fill_value=0,
    ).reset_index()

    # make source cols int
    for col in ["af", "ct_ai", "ct_noai", "other"]:
        authors[col] = authors[col].astype(int)

    authors = authors.merge(strength_labels, on=["author", "quarter"], how="left")

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
    slices = [filter_batches[i : i + 50] for i in range(0, len(filter_batches), 50)]

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

    institutions = _institution_preprocessing(institutions)

    # process pdb data
    author_submissions = _process_pdb_data(pdb_submissions)

    # drop publication_date from pdb_submissions
    pdb_submissions = pdb_submissions.drop(columns=["publication_date"])

    mesh_terms["term_group"] = mesh_terms["tree_number"].apply(
        lambda x: str(x)[:1] if x is not None else None
    )
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
        data = define_high_pdb_authors(data, author_submissions)

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

        for col in ["concepts", "mesh_terms", "grants", "topics", "ids"]:
            try:
                data.drop(
                    columns=[col],
                    inplace=True,
                )
            except KeyError:
                pass

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
        cols_to_fill = [
            col for col in candidate_authors.columns if col not in ["author", "quarter"]
        ]
        for col in cols_to_fill:
            data[col] = data.groupby("author")[col].ffill().fillna(0).astype(int)

        # get patent citations
        data = get_patent_citations(data, patents_data)

        # get pdb metrics
        data = data.merge(
            pdb_submissions, on="id", how="left", indicator=True
        )  # TO DOUBLECHECK
        data["pdb_submission"] = data["_merge"].apply(
            lambda x: True if x == "both" else False
        )
        data = data.drop(columns=["_merge"])

        # get icite_counts
        icite_outputs = create_cc_counts(data[["id", "doi", "pmid"]], icite_data)
        data = data.merge(icite_outputs, on="id", how="left")

        # concatenate the data
        output_data = pd.concat([output_data, data], ignore_index=True)

    # create top quantile pre2021 pdb activity
    # get one row per author with their pre2021 pdb submissions
    author_pdb = output_data[
        ["author", "num_pdb_submissions_pre2021"]
    ].drop_duplicates()
    # calculate threshold based on author-level data
    pdb_threshold = author_pdb["num_pdb_submissions_pre2021"].quantile(0.75)
    # merge back to full dataset
    output_data["high_pdb_pre2021"] = (
        output_data["num_pdb_submissions_pre2021"] > pdb_threshold
    )

    return output_data


def _institution_preprocessing(institutions: pd.DataFrame) -> pd.DataFrame:
    institutions = institutions.drop(columns=["author"]).drop_duplicates(
        subset=["institution"]
    )
    institutions.rename(
        columns={"cited_by_count": "institution_cited_by_count"}, inplace=True
    )

    return institutions


def define_high_pdb_authors(
    data: pd.DataFrame, author_submissions: pd.DataFrame
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

    author_submissions_pre2021 = author_submissions[
        author_submissions["publication_date"] < "2021-01-01"
    ]
    data_db = (
        data[["author"]]
        .drop_duplicates(subset=["author"])
        .merge(author_submissions_pre2021[["author", "id"]], on="author", how="left")
        .drop_duplicates(subset=["author", "id"])
    )

    # create a column that is 1 if notna else 0
    data_db["num_pdb_submissions_pre2021"] = data_db["id"].notna().astype(int)

    # group by author and sum pdb_submission
    data_db = (
        data_db.groupby("author")["num_pdb_submissions_pre2021"].sum().reset_index()
    )

    # Merge data_db to include pdb_submission column
    data = data.merge(
        data_db[["author", "num_pdb_submissions_pre2021"]],
        on="author",
        how="left",
    )

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
        lambda x: [y.get("descriptor_ui") for y in x] if x is not None else []
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

    # fill NaN for columns that are mesh_
    mesh_cols = [col for col in data.columns if col.startswith("mesh_")]
    data[mesh_cols] = data[mesh_cols].fillna(0)

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

    data_2020 = data[
        (data["publication_date"] >= "2020-01-01")
        & (data["publication_date"] <= "2020-12-31")
    ]

    # explode concepts
    data_2020 = data_2020.explode("concepts")

    # get the first element of each list
    data_2020["concepts"] = data_2020["concepts"].apply(
        lambda x: x[0] if isinstance(x, np.ndarray) and len(x) > 0 else None
    )

    data_2020["covid_2020"] = data_2020["concepts"].apply(lambda x: (x == "C524204448"))

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

    # fill NaN with 0
    data["covid_share_2020"] = data["covid_share_2020"].fillna(0)

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

    def safe_mode(series):
        mode = series.mode()
        return mode.iloc[0] if not mode.empty else np.nan

    # aggregation dictionary
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
        "num_uniprot_structures": ("num_uniprot_structures", "sum"),
        "num_pdb_ids": ("num_pdb_ids", "sum"),
        "num_primary_submissions": ("num_primary_submissions", "sum"),
        "score_mean": ("score_mean", "mean"),
        "complexity_sum": ("complexity_sum", "sum"),
        "complexity_mean": ("complexity_mean", "mean"),
        "organism_rarity_mean": ("organism_rarity_mean", "mean"),
        "organism_rarity_max": ("organism_rarity_max", "mean"),
        "num_diseases": ("num_diseases", "sum"),
        "resolution_mean": ("resolution_mean", "mean"),
        "R_free_mean": ("R_free_mean", "mean"),
        "pdb_submission": ("pdb_submission", "sum"),
        "mean_tmscore": ("mean_tmscore", "mean"),
        "max_tmscore": ("max_tmscore", "max"),
        "normalised_mean_tmscore": ("normalised_mean_tmscore", "mean"),
        "normalised_max_tmscore": ("normalised_max_tmscore", "max"),
        "institution": ("institution", "first"),
        "institution_cited_by_count": ("institution_cited_by_count", "first"),
        "institution_country_code": ("country_code", "first"),
        "institution_type": ("type", "first"),
        "institution_2yr_mean_citedness": ("2yr_mean_citedness", "first"),
        "institution_h_index": ("h_index", "first"),
        "institution_i10_index": ("i10_index", "first"),
        "institution_works_count": ("works_count", "first"),
        "af": ("af", safe_mode),
        "ct_ai": ("ct_ai", safe_mode),
        "ct_noai": ("ct_noai", safe_mode),
        "other": ("other", safe_mode),
        "af_mixed": ("af_mixed", safe_mode),
        "af_strong": ("af_strong", safe_mode),
        "af_unknown": ("af_unknown", safe_mode),
        "af_weak": ("af_weak", safe_mode),
        "ct_ai_mixed": ("ct_ai_mixed", safe_mode),
        "ct_ai_strong": ("ct_ai_strong", safe_mode),
        "ct_ai_unknown": ("ct_ai_unknown", safe_mode),
        "ct_ai_weak": ("ct_ai_weak", safe_mode),
        "ct_noai_mixed": ("ct_noai_mixed", safe_mode),
        "ct_noai_strong": ("ct_noai_strong", safe_mode),
        "ct_noai_unknown": ("ct_noai_unknown", safe_mode),
        "ct_noai_weak": ("ct_noai_weak", safe_mode),
        "other_mixed": ("other_mixed", safe_mode),
        "other_strong": ("other_strong", safe_mode),
        "other_unknown": ("other_unknown", safe_mode),
        "other_weak": ("other_weak", safe_mode),
        "af_with_intent": ("af_with_intent", safe_mode),
        "ct_ai_with_intent": ("ct_ai_with_intent", safe_mode),
        "ct_noai_with_intent": ("ct_noai_with_intent", safe_mode),
        "other_with_intent": ("other_with_intent", safe_mode),
        "primary_field": ("primary_field", safe_mode),
        "author_position": ("author_position", safe_mode),
        "high_pdb_pre2021": ("high_pdb_pre2021", safe_mode),
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


def create_cc_counts(data, icite_data):
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
    icite_data["pmid"] = icite_data["pmid"].astype(str)
    icite_data["doi"] = icite_data["doi"].astype(str)

    # merge on 'pmid'
    data_pmid = data.merge(icite_data[["pmid", "cited_by_clin"]], how="left", on="pmid")
    data_pmid = data_pmid.drop_duplicates(subset=["id"])

    # merge on 'doi'
    data_doi = data.merge(icite_data[["doi", "cited_by_clin"]], how="left", on="doi")
    data_doi = data_doi.drop_duplicates(subset=["id"])

    combined_data = pd.concat([data_pmid, data_doi]).drop_duplicates(subset=["id"])

    logger.info("Exploding cited_by_clin")
    combined_data = combined_data[combined_data["cited_by_clin"].astype(str) != "nan"]

    # drop if cited_by_clin is None
    combined_data = combined_data[combined_data["cited_by_clin"].notnull()]

    combined_data["cited_by_clin"] = combined_data["cited_by_clin"].apply(
        lambda x: x.split(" ")
    )

    # create count column
    combined_data["ca_count"] = combined_data["cited_by_clin"].apply(len)

    # merge back with data, fill with 0 for missing count
    data = data.merge(combined_data[["id", "ca_count"]], on="id", how="left")
    data["ca_count"] = data["ca_count"].fillna(0).astype(int)

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
    data["pmid"] = data["ids"].apply(
        lambda x: (
            x.get("pmid", None).replace("https://pubmed.ncbi.nlm.nih.gov/", "")
            if x.get("pmid")
            else None
        )
    )
    data["pmcid"] = data["ids"].apply(
        lambda x: (
            x.get("pmcid", None).replace(
                "https://www.ncbi.nlm.nih.gov/pmc/articles/", ""
            )
            if x.get("pmcid")
            else None
        )
    )
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


def _process_pdb_data(pdb_submissions: pd.DataFrame) -> pd.DataFrame:
    return (
        pdb_submissions.explode("authorships")
        .assign(
            authorships=lambda x: x["authorships"].apply(
                lambda y: y[0] if isinstance(y, np.ndarray) and len(y) > 0 else None
            )
        )
        .rename(columns={"authorships": "author"})
    )


def get_patent_citations(
    data: pd.DataFrame, patents_data: pd.DataFrame
) -> pd.DataFrame:
    """
    Get patent citations for the given data using multiple identifier types (DOI, PMID, PMCID).

    Args:
        data (pandas.DataFrame): The input data containing paper identifiers.
        patents_data (pandas.DataFrame): The patent data with NPL citations.

    Returns:
        pandas.DataFrame: The data with patent citations merged and aggregated.
    """

    # remove prefix
    data["doi"] = data["doi"].replace("https://doi.org/", "", regex=True)

    # Create separate matches for each ID type
    matches = []

    # Match DOIs
    doi_matches = data.merge(
        patents_data, left_on="doi", right_on="NPL Resolved External ID(s)", how="inner"
    )
    matches.append(doi_matches)

    # Match PMIDs
    pmid_matches = data.merge(
        patents_data,
        left_on="pmid",
        right_on="NPL Resolved External ID(s)",
        how="inner",
    )
    matches.append(pmid_matches)

    # Match PMCIDs
    pmcid_matches = data.merge(
        patents_data,
        left_on="pmcid",
        right_on="NPL Resolved External ID(s)",
        how="inner",
    )
    matches.append(pmcid_matches)

    # Combine all matches and remove duplicates
    patent_matches = pd.concat(matches, ignore_index=True)
    patent_matches = patent_matches.drop_duplicates(subset=["id", "Title"])

    # Group by paper ID and aggregate patent information
    patent_matches_grouped = (
        patent_matches.groupby(["id"])
        .agg(
            patent_count=pd.NamedAgg(column="Title", aggfunc="count"),
            CPCs=pd.NamedAgg(
                column="CPC Classifications",
                aggfunc=lambda x: ";;".join(map(str, x)),
            ),
            patent_citation=pd.NamedAgg(column="Cited by Patent Count", aggfunc="sum"),
        )
        .reset_index()
    )

    # Merge with original data and fill missing values
    data = data.merge(patent_matches_grouped, on="id", how="left")
    data["patent_count"] = data["patent_count"].fillna(0)

    return data
