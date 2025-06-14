"""
This module contains utility functions for the data collection pipeline.
"""

import logging
from typing import List, Tuple
import requests
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def get_sb_candidate_authors(data: pd.DataFrame) -> List[Tuple[str, str]]:
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

    # instead return all author, institution pairs when position is last
    candidate_authors = (
        author_data[author_data["position"] == "last"]
        .groupby(["author", "institution"])["id"]
        .nunique()
    )

    return list(candidate_authors.index)


def separate_ct_from_seed(
    baseline_data,
    seed_baseline_data,
    ct_data,
):
    """Separate the CT data from the seed data into three technology types: ai, pp, and sb."""
    # create ct_level_0_ai_papers
    ct_level_0_ai_papers = baseline_data[
        baseline_data["parent_id"].isin(
            ct_data[ct_data["label"].str.contains("ai")]["parent_id"]
        )
    ]

    # create ct_level_0_pp_papers
    ct_level_0_pp_papers = baseline_data[
        baseline_data["parent_id"].isin(
            ct_data[ct_data["label"].str.contains("pp")]["parent_id"]
        )
    ]

    # create ct_level_0_sb_papers
    ct_level_0_sb_papers = baseline_data[
        baseline_data["parent_id"].isin(
            ct_data[ct_data["label"].str.contains("sb")]["parent_id"]
        )
    ]

    # subset the seed data for rows that do not have the id in ct_data (as parent_id)
    other_level_0_papers = seed_baseline_data[
        ~seed_baseline_data["id"].isin(ct_data["parent_id"])
    ]

    return (
        ct_level_0_ai_papers,
        ct_level_0_pp_papers,
        ct_level_0_sb_papers,
        other_level_0_papers,
    )


def compute_apf(author_data: pd.DataFrame) -> pd.DataFrame:
    """
    Compute the Average Position Factor (APF) for each author, institution pair
        based on the given author data.

    Parameters:
        author_data (pd.DataFrame): A DataFrame containing author data with
            columns 'author', 'institution', 'position', and 'id'.

    Returns:
        pd.DataFrame: A DataFrame containing the APF values for each author,
            institution pair, along with the share of unique id they appear as
            first, middle, or last author.
    """

    # for each author, institution unique pair,
    total_counts = author_data.groupby(["author", "institution", "year"])[
        "id"
    ].nunique()
    position_counts = (
        author_data.groupby(["author", "institution", "year", "position"])["id"]
        .nunique()
        .unstack()
    )
    position_shares = position_counts.div(total_counts, axis=0)

    # create the three columns, fill nan with 0
    position_shares = position_shares.fillna(0).reset_index()
    position_shares.columns = [
        "author",
        "institution",
        "year",
        "first",
        "middle",
        "last",
    ]

    # compute the average position factor (APF) for each author, institution pair.
    position_shares["apf"] = (
        position_shares["last"] + 0.5 * position_shares["first"]
    ) * (1 - position_shares["middle"])

    return position_shares


def compute_avg_citation_count(author_data: pd.DataFrame) -> pd.DataFrame:
    """
    Compute the average citation count for each author, institution pair
        based on the given author data.

    Parameters:
        author_data (pd.DataFrame): A DataFrame containing author data with
            columns 'author', 'institution', 'position', and 'id'.

    Returns:
        pd.DataFrame: A DataFrame containing the average citation count for each author,
            institution pair.
    """

    # for each author, institution unique pair,
    avg_citation_count = (
        author_data.groupby(["author", "institution", "year"])["cited_by_count"]
        .mean()
        .apply(np.log1p)
        .reset_index()
    )

    return avg_citation_count


def compute_sample_publication_count(author_data: pd.DataFrame) -> pd.DataFrame:
    """
    Compute the sample publication count for each author, institution pair
        based on the given author data.

    Parameters:
        author_data (pd.DataFrame): A DataFrame containing author data with
            columns 'author', 'institution', 'position', and 'id'.

    Returns:
        pd.DataFrame: A DataFrame containing the sample publication count for each author,
            institution pair.
    """

    # for each author, institution unique pair,
    sample_publication_count = (
        author_data.groupby(["author", "institution", "year"])["id"]
        .nunique()
        .apply(np.log1p)
        .reset_index()
    )

    return sample_publication_count


def appears_in_three_consecutive_years(group: pd.DataFrame) -> bool:
    """
    Checks if a group of data appears in three consecutive years.

    Args:
        group (pandas.DataFrame): The group of data to check.

    Returns:
        bool: True if the group appears in three consecutive years, False otherwise.
    """
    years = group["year"].sort_values().diff().eq(1)
    return years.rolling(3).sum().eq(3).any()


def is_likely_pi(group: pd.DataFrame, quantile: float) -> bool:
    """
    Determines if a group is likely to be a principal investigator (PI) based
        on the given criteria.

    Args:
        group (pandas.DataFrame): A DataFrame representing a group of data.
        quantile (float): The quantile threshold for the 'score' column.

    Returns:
        bool: True if the group is likely to be a PI, False otherwise.
    """
    group = group.sort_values("year")
    above = group["score"] > quantile
    return (
        above.rolling(3).sum().eq(3).any()
        or (group["year"].isin([2022, 2023]) & above).any()
    )


def fetch_institution_data(institution):
    """
    Fetches data for a given institution from the OpenAlex API.

    Args:
        institution (str): The name of the institution.

    Returns:
        dict: A dictionary containing the fetched institution data, including the institution name,
            country code, type, works count, cited by count, 2-year mean citedness, h-index, and
            i10-index.

        If the request fails or the institution data is not found, None is returned.
    """
    url = f"https://api.openalex.org/institutions/{institution}"
    while True:
        try:
            response = requests.get(url, timeout=40)

            # Process the response if the request was successful
            institution_data = response.json()

            # extract country_code, type, works_count, cited_by_count, summary_stats
            institution_data = {
                "institution": institution,
                "country_code": institution_data.get("country_code"),
                "type": institution_data.get("type"),
                "works_count": institution_data.get("works_count"),
                "cited_by_count": institution_data.get("cited_by_count"),
                "2yr_mean_citedness": institution_data.get("summary_stats", {}).get(
                    "2yr_mean_citedness"
                ),
                "h_index": institution_data.get("summary_stats", {}).get("h_index"),
                "i10_index": institution_data.get("summary_stats", {}).get("i10_index"),
            }
            break
        except:  # pylint: disable=bare-except
            logger.info("Error fetching for institution %s", institution)
            return None
    return institution_data


def explode_author_data(data: pd.DataFrame) -> pd.DataFrame:
    """
    Explodes the author data in the given DataFrame and returns a new DataFrame
    with the exploded author data.

    Args:
        data (pd.DataFrame): The input DataFrame containing the author data.

    Returns:
        pd.DataFrame: The new DataFrame with exploded author data, including columns for
        id, author, institution, and position.
    """

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

    # drop author A9999999999 (placeholder for missing author data)
    author_data = author_data[author_data["author"] != "A9999999999"]

    # drop "" institution (placeholder for missing institution data)
    author_data = author_data[author_data["institution"] != ""]

    # get author year
    author_data["year"] = author_data["publication_date"].str[:4]

    return author_data[
        ["id", "year", "author", "institution", "position", "cited_by_count"]
    ]
