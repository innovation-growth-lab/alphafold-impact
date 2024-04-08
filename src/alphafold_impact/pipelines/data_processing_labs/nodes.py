"""
This is a boilerplate pipeline 'data_processing_labs'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict, Generator
import pandas as pd
from kedro.io import AbstractDataset

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
