"""
This is a boilerplate pipeline 'data_processing_labs'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict, Generator
import pandas as pd
import numpy as np
from kedro.io import AbstractDataset
from sklearn.preprocessing import MinMaxScaler

logger = logging.getLogger(__name__)


def _explode_author_data(data: pd.DataFrame) -> pd.DataFrame:
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


def _compute_apf(author_data: pd.DataFrame) -> pd.DataFrame:
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


def _compute_avg_citation_count(author_data: pd.DataFrame) -> pd.DataFrame:
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
        .reset_index()
    )

    return avg_citation_count


def _compute_sample_publication_count(author_data: pd.DataFrame) -> pd.DataFrame:
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
        .reset_index()
    )

    return sample_publication_count


def calculate_lab_determinants(
    dict_loader: AbstractDataset,
) -> Generator[Dict[str, pd.DataFrame], None, None]:
    """
    This node calculates the lab determinants for each lab based on the given data.

    Parameters:
        dict_loader (AbstractDataset): The input data loader containing the lab data.

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

        # break atuhorship nested dictionary jsons, create triplets of authorship
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

        # for each dictionary, only preserve
        author_data = _explode_author_data(data)

        # lets drop authors whose row count is less than 10 (reduces chances of reassessment)
        author_data = author_data.groupby(["author", "institution"]).filter(
            lambda x: len(x) >= 10
        )

        # compute average annual citation count
        avg_citation_count = _compute_avg_citation_count(author_data)

        # compute sample publication count
        sample_publication_count = _compute_sample_publication_count(author_data)

        # compute apf
        apf_data = _compute_apf(author_data)

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


def assign_lab_label(
    candidate_data: pd.DataFrame,
    ground_truth_data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Assigns lab labels to candidate data based on matching with ground truth data.

    Args:
        candidate_data (pd.DataFrame): The candidate data containing author information.
        ground_truth_data (pd.DataFrame): The ground truth data containing PI IDs.

    Returns:
        pd.DataFrame: The candidate data with lab labels assigned.

    """
    logger.info("Cleaning NIH PI data")
    ground_truth_data["pi_id"] = ground_truth_data["pi_id"].str.replace(
        "https://openalex.org/", ""
    )

    logger.info("Match author in candidate data")
    candidate_data["label"] = (
        candidate_data["author"].isin(ground_truth_data["pi_id"]).astype(int)
    )
    # we only match 10 unique researchers

    logger.info("Sorting candidate data")
    candidate_data = candidate_data.sort_values(by=["author", "year"]).reset_index(
        drop=True
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

    logger.info("Calculating quantiles")
    quantiles = candidate_data["mean_apf"].quantile([0.25, 0.5, 0.75])

    logger.info(
        "Filtering candidate data based on having at least one year above 75th percentile of APF"
    )
    candidate_data = candidate_data[
        candidate_data.groupby(["author", "institution"])["apf"].transform(
            lambda x: x.max() > quantiles[0.75]
        )
    ]

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
    candidate_data[["cited_by_count", "publication_count"]] = scaler.fit_transform(
        candidate_data[["cited_by_count", "publication_count"]]
    )

    # assign weights
    weights = {
        "apf": 0.4,
        "mean_apf": 0.1,
        "apf_3yr_avg": 0.2,
        "publication_count": 0.2,
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
    quantile_90 = candidate_data["score"].quantile(0.90)
    likely_pis = candidate_data.groupby(["author", "institution"]).filter(
        is_likely_pi, quantile=quantile_90
    )

    return likely_pis
