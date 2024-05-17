"""
This is a boilerplate pipeline 'data_analysis_ecr'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
from typing import Generator
from kedro.io import AbstractDataset

logger = logging.getLogger(__name__)


def run_analysis(data_loaders: AbstractDataset, af_data: pd.DataFrame) -> Generator:
    """
    Runs the analysis on the provided data loaders and AlphaFold data.

    Args:
        data_loaders (AbstractDataset): The data loaders to iterate over.
        af_data (pd.DataFrame): The AlphaFold data.

    Yields:
        dict: A dictionary containing the computed analysis results.

    """
    af_data = af_data[af_data["level"] == "0"]

    af_sf = []
    no_af_sf = []
    af_f = []
    no_af_f = []
    af_counts = []
    no_counts = []
    af_citations = []
    no_citations = []

    for i, loader in enumerate(data_loaders.values()):
        logger.info("Loading data batch %d / %d", i + 1, len(data_loaders))
        dta = loader()
        dta = dta[
            [
                "id",
                "doi",
                "publication_date",
                "authorships",
                "concepts",
                "topics",
                "grants",
                "t0",
                "t1",
                "t2",
                "t3",
                "t4",
                "t5",
                "t6",
            ]
        ]

        # preprocess data
        dta = _preprocess(dta)

        logger.info("Creating data with and without AlphaFold data")
        authors_with_af, authors_with_no_af = _create_af_and_no_af_data(dta, af_data)

        logger.info("Computing subfield shares")
        af_sf_dta, no_sf_dta = _compute_topic_shares(
            authors_with_af, authors_with_no_af, "subfields"
        )

        logger.info("Computing field shares")
        af_f_dta, no_f_dta = _compute_topic_shares(
            authors_with_af, authors_with_no_af, "fields"
        )

        logger.info("Computing publication counts")
        af_counts_dta = _compute_publication_count(authors_with_af)
        no_counts_dta = _compute_publication_count(authors_with_no_af)

        logger.info("Computing average citations")
        af_citations_dta = _compute_average_citations(authors_with_af)
        no_citations_dta = _compute_average_citations(authors_with_no_af)

        af_sf.append(af_sf_dta)
        no_af_sf.append(no_sf_dta)
        af_f.append(af_f_dta)
        no_af_f.append(no_f_dta)
        af_counts.append(af_counts_dta)
        no_counts.append(no_counts_dta)
        af_citations.append(af_citations_dta)
        no_citations.append(no_citations_dta)

    af_sf = pd.concat(af_sf, ignore_index=True)
    no_af_sf = pd.concat(no_af_sf, ignore_index=True)
    af_f = pd.concat(af_f, ignore_index=True)
    no_af_f = pd.concat(no_af_f, ignore_index=True)
    af_counts = pd.concat(af_counts, ignore_index=True)
    no_counts = pd.concat(no_counts, ignore_index=True)
    af_citations = pd.concat(af_citations, ignore_index=True)
    no_citations = pd.concat(no_citations, ignore_index=True)

    return (
        af_sf,
        no_af_sf,
        af_f,
        no_af_f,
        af_counts,
        no_counts,
        af_citations,
        no_citations,
    )


def _preprocess(data: pd.DataFrame):
    """
    Preprocesses the given DataFrame by performing the following operations:
    - Renames the 'authorships' column to 'author'
    - Removes the 'https://openalex.org/' prefix from the 'id' column values
    - Extracts the subfields from the 'topics' column and stores them in a new 'subfields' column

    Args:
        data (pd.DataFrame): The input DataFrame to be preprocessed.

    Returns:
        pd.DataFrame: The preprocessed DataFrame.
    """
    data.rename(columns={"authorships": "author"}, inplace=True)
    data["id"] = data["id"].str.replace("https://openalex.org/", "")
    data["subfields"] = data["topics"].apply(
        lambda x: [y[3] for y in x] if x is not None and len(x) > 0 else []
    )
    data["fields"] = data["topics"].apply(
        lambda x: [y[5] for y in x] if x is not None and len(x) > 0 else []
    )
    return data


def _create_af_and_no_af_data(data: pd.DataFrame, af_data: pd.DataFrame):
    af_matched_data = data[["id", "author"]].merge(af_data, on="id", how="inner")
    af_matched_data = af_matched_data.sort_values("publication_date").drop_duplicates(
        "author"
    )
    author_to_pub_date_af = (
        af_matched_data[["author", "publication_date"]]
        .set_index("author")
        .to_dict()["publication_date"]
    )

    # create a dataframe with rows that have af_pub_date
    authors_with_af = data[data["author"].isin(author_to_pub_date_af.keys())]
    authors_with_no_af = data[~data["author"].isin(author_to_pub_date_af.keys())]

    authors_with_af["af_pub_date"] = authors_with_af["author"].map(
        author_to_pub_date_af
    )
    authors_with_no_af["af_pub_date"] = "2022-01-01"

    return authors_with_af, authors_with_no_af


def _compute_topic_shares(authors_with_af, authors_with_no_af, column):
    """
    Computes the subfield shares for authors with and without AlphaFold data.

    Args:
        data (pd.DataFrame): The input DataFrame containing author data.
        af_data (pd.DataFrame): The DataFrame containing AlphaFold data.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames.
            The first DataFrame represents the subfield shares for authors with AlphaFold data.
            The second DataFrame represents the subfield shares for authors without AlphaFold data.
    """

    af_using_shares = _topic_shares(authors_with_af, column=column)
    no_af_using_shares = _topic_shares(authors_with_no_af, column=column)

    return af_using_shares, no_af_using_shares


def _topic_shares(df, column: str):
    """
    Calculate the share of each subfield for pre- and post-AlphaFold data.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data.

    Returns:
        pandas.DataFrame: The resulting DataFrame with the calculated shares.

    """
    pre_af_data = df[df["publication_date"] < df["af_pub_date"]]
    post_af_data = df[df["publication_date"] >= df["af_pub_date"]]

    # calculate the share of each subfield for pre-AF and post-AF data
    pre_af_share = _calculate_topic_share(pre_af_data, column)
    post_af_share = _calculate_topic_share(post_af_data, column)

    # join the two DataFrames together
    pre_af_share["status"] = "pre"
    post_af_share["status"] = "post"
    result = pd.concat([pre_af_share, post_af_share])
    result.fillna(0, inplace=True)

    # sort columns, status first and then alphabetically
    result = result[["status"] + sorted(result.columns.difference(["status"]))]

    return result


def _calculate_topic_share(df, column: str):
    """
    Calculate the share of each subfield for each author in a DataFrame.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data.

    Returns:
        pandas.DataFrame: The resulting DataFrame with the share of each subfield for each author.
    """

    # explode the subfields column into multiple rows
    df = df.explode(column)
    # group by author and subfields and calculate the count
    df = df.groupby(["author", column]).size().reset_index(name="count")
    # calculate the total count for each author
    total_count = df.groupby("author")["count"].sum()
    # calculate the share of each subfield
    df["share"] = df.apply(
        lambda row: row["count"] / total_count[row["author"]], axis=1
    )
    # pivot the DataFrame to get one column for each subfield
    df = df.pivot(index="author", columns=column, values="share")
    return df


def _compute_publication_count(df):
    """
    Compute the publication count for each author, categorized as pre or post AlphaFold.

    Args:
        df (pandas.DataFrame): The input DataFrame containing publication data.

    Returns:
        pandas.DataFrame: A DataFrame with the publication count for each author, categorized as
            pre or post AlphaFold.
                The DataFrame has columns 'author', 'status', and 'count'.
    """

    # add a column 'status' that indicates whether each publication is pre or post AlphaFold
    df["status"] = df.apply(
        lambda row: "pre" if row["publication_date"] < row["af_pub_date"] else "post",
        axis=1,
    )

    # group by author and status, and count the number of publications
    publication_count = (
        df.groupby(["author", "status"]).size().reset_index(name="count")
    )

    # create a MultiIndex with all combinations of authors and periods
    all_authors = df["author"].unique()
    all_periods = ["pre", "post"]
    index = pd.MultiIndex.from_product(
        [all_authors, all_periods], names=["author", "status"]
    )

    # reindex the DataFrame with the new index and fill NaN values with 0
    publication_count = (
        publication_count.set_index(["author", "status"])
        .reindex(index, fill_value=0)
        .reset_index()
    )

    # pivot the DataFrame to get one column for pre and post AlphaFold
    publication_count = publication_count.pivot(
        index="author", columns="status", values="count"
    )

    # fill NaN values with 0
    publication_count.fillna(0, inplace=True)

    # reshape the DataFrame from wide format to long format
    publication_count = publication_count.reset_index().melt(
        id_vars="author", var_name="status", value_name="count"
    )

    return publication_count


def _compute_average_citations(df):
    """
    Compute the average number of t0 and t1 citations for each author, categorized as pre
    or post AlphaFold.

    Args:
        df (pandas.DataFrame): The input DataFrame containing publication data.

    Returns:
        pandas.DataFrame: A DataFrame with the average number of t0 and t1 citations for each
         author, categorized as pre or post AlphaFold. The DataFrame has columns 'author', 'status',
         't0_citations', and 't1_citations'.
    """

    # add a column 'status' that indicates whether each publication is pre or post AlphaFold
    df["status"] = df.apply(
        lambda row: "pre" if row["publication_date"] < row["af_pub_date"] else "post",
        axis=1,
    )

    # group by author and status, and calculate the mean number of t0 and t1 citations
    average_citations = (
        df.groupby(["author", "status"])[["t0", "t1"]].mean().reset_index()
    )

    # create a MultiIndex with all combinations of authors and periods
    all_authors = df["author"].unique()
    all_periods = ["pre", "post"]
    index = pd.MultiIndex.from_product(
        [all_authors, all_periods], names=["author", "status"]
    )

    # reindex the DataFrame with the new index and fill NaN values with 0
    average_citations = (
        average_citations.set_index(["author", "status"])
        .reindex(index, fill_value=0)
        .reset_index()
    )

    # pivot the DataFrame to get one column for pre and post AlphaFold for each citation type
    average_citations = average_citations.melt(
        id_vars=["author", "status"],
        var_name="citation_type",
        value_name="average_citations",
    )

    return average_citations
