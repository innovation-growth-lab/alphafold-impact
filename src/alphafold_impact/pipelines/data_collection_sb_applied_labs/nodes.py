"""
This is a boilerplate pipeline 'data_collection_other_labs'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict, List, Tuple, Generator
from kedro.io import AbstractDataset
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

from ..data_processing_labs.nodes import (  # pylint: disable=E0402
    _explode_author_data,
    _compute_avg_citation_count,
    _compute_sample_publication_count,
    _compute_apf,
    appears_in_three_consecutive_years,
    is_likely_pi,
)

logger = logging.getLogger(__name__)


def _get_applied_candidate_authors(data: pd.DataFrame) -> List[Tuple[str, str]]:
    """
    Get the last authors of the given papers.

    Args:
        data (pd.DataFrame): The input DataFrame containing the papers.

    Returns:
        List[Tuple[str, str]]: The list of last authors.
    """
    # filter data to only consider levels 1 and 2
    applied_ids = data[data["level"].astype(str).isin(["1", "2"])]['id'].unique()

    # Identify the 'id' values in 'data' where 'level' is in the other levels ("0", "seed")
    sb_ids = data[~data["level"].astype(str).isin(["1", "2"])]['id'].unique()

    # Select the rows in 'applied_ids' that are not in 'other_ids'
    data = data[data['id'].isin(applied_ids) & ~data['id'].isin(sb_ids)]

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

    # calculate counts for each author, institution pair
    counts = (
        author_data[author_data["position"] == "last"]
        .groupby(["author", "institution"])["id"]
        .count()
    )

    # determine the quantile threshold
    threshold = counts.quantile(0.75)

    # filter pairs based on the threshold
    alphafold_authors = counts[counts >= threshold]

    return list(alphafold_authors.index)


def get_candidate_authors(
    alphafold_data: pd.DataFrame,
    ct_data: pd.DataFrame,
    other_data: pd.DataFrame,
):
    """
    Get a list of candidate authors from different data sources.

    Args:
        alphafold_data (pd.DataFrame): Dataframe containing AlphaFold data.
        ct_data (pd.DataFrame): Dataframe containing CT data.
        other_data (pd.DataFrame): Dataframe containing other data.

    Returns:
        list: A list of candidate authors from different data sources.
    """
    # get the last authors from the AlphaFold data
    alphafold_authors = _get_applied_candidate_authors(alphafold_data)

    # get the last authors from the CT data
    ct_authors = _get_applied_candidate_authors(ct_data)

    # get the last authors from the other SB data
    other_authors = _get_applied_candidate_authors(other_data)

    # create a unique list of authors
    authors = list(set(alphafold_authors + ct_authors + other_authors))

    return authors, alphafold_authors, ct_authors, other_authors


def create_candidates_map(
    alphafold_authors: List,
    ct_authors: List,
    other_authors: List,
) -> pd.DataFrame:
    """
    Create a candidates map by combining the author data from different sources.

    Args:
        alphafold_authors (List): A list of authors from the AlphaFold source.
        ct_authors (List): A list of authors from the CT source.
        other_authors (List): A list of authors from other sources.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the combined author data.

    """

    # create pandas dataframes for alphafold_authors, ct_authors, other_authors
    alphafold_df = pd.DataFrame(alphafold_authors, columns=["author", "institution"])
    alphafold_df["seed"] = "alphafold"

    ct_df = pd.DataFrame(ct_authors, columns=["author", "institution"])
    ct_df["seed"] = "ct"

    other_df = pd.DataFrame(other_authors, columns=["author", "institution"])
    other_df["seed"] = "other"

    # combine the candidate dataframes
    candidate_data = pd.concat([alphafold_df, ct_df, other_df])

    # check for duplicates in candidate_data
    # candidate_data = candidate_data.drop_duplicates(subset=["author", "institution"])

    candidate_data["seed"] = candidate_data.groupby(["author", "institution"])[
        "seed"
    ].transform(lambda x: ",".join(x.unique()))
    candidate_data = candidate_data.drop_duplicates(subset=["author", "institution"])

    # ad-hoc fixes
    # if "seed" is alphafold,other let's go back to alphafold
    candidate_data["seed"] = candidate_data["seed"].apply(
        lambda x: "alphafold" if x == "alphafold,other" else x
    )

    # if "seed" is alphafold,ct,other let's go to alphafold,ct
    candidate_data["seed"] = candidate_data["seed"].apply(
        lambda x: "alphafold,ct" if x == "alphafold,ct,other" else x
    )

    # if ct,other let's go back to ct
    candidate_data["seed"] = candidate_data["seed"].apply(
        lambda x: "ct" if x == "ct,other" else x
    )

    return candidate_data


def calculate_lab_determinants(
    dict_loader: AbstractDataset,
    alphafold_authors: List,
    ct_authors: List,
    other_authors: List,
) -> Generator[Dict[str, pd.DataFrame], None, None]:
    """
    This node calculates the lab determinants for each lab based on the given data.

    Parameters:
        dict_loader (AbstractDataset): The input data loader containing the lab data.
        alphafold_authors (List): The list of AlphaFold authors.
        ct_authors (List): The list of CT authors.
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

        # get mapping between candidates and chain seeds
        candidate_map = create_candidates_map(
            alphafold_authors, ct_authors, other_authors
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


def assign_lab_label(
    candidate_data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Assigns lab labels to candidate data based on matching with ground truth data.

    Args:
        candidate_data (pd.DataFrame): The candidate data containing author information.

    Returns:
        pd.DataFrame: The candidate data with lab labels assigned.

    """

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
    quantile_75 = candidate_data["score"].quantile(0.75)
    likely_pis = candidate_data.groupby(["author", "institution"]).filter(
        is_likely_pi, quantile=quantile_75
    )

    return likely_pis
