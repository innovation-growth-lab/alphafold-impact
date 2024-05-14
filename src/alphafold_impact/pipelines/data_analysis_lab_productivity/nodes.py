"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
import numpy as np
from kedro.io import AbstractDataset

logger = logging.getLogger(__name__)


def get_sb_lab_outputs(
    data_loaders: AbstractDataset,
    mapping_df: AbstractDataset,
):
    """
    Retrieves outputs from data loaders and performs data processing.

    Args:
        data_loaders (AbstractDataset): A dictionary-like object containing data loaders.
        mapping_df (AbstractDataset): A dataset containing mapping information.

    Returns:
        pd.DataFrame: The processed data.

    """
    outputs = []
    for i, loader in enumerate(data_loaders.values()):
        logger.info("Loading data batch %d / %d", i + 1, len(data_loaders))
        data_batch = loader()

        # drop "display_name", "authorships"
        data_batch = data_batch.drop(
            columns=["display_name", "authorships", "counts_by_year", "ids"]
        )

        # transform mesh, concepts, topics to be a list of the ids (ie first item in sublists)
        data_batch["mesh_terms"] = data_batch["mesh_terms"].apply(
            lambda x: [y[0] for y in x] if x is not None else []
        )
        data_batch["concepts"] = data_batch["concepts"].apply(
            lambda x: [y[0] for y in x] if x is not None else []
        )
        data_batch["topics"] = data_batch["topics"].apply(
            lambda x: [y[0] for y in x] if x is not None else []
        )

        outputs.append(data_batch)
    data = pd.concat(outputs)

    logger.info("Merging data with mapping information")
    mapping_df.drop_duplicates(subset="author", inplace=True)
    data = data.merge(
        mapping_df[["author", "seed"]], left_on="pi_id", right_on="author", how="left"
    )

    data["publication_date"] = pd.to_datetime(data["publication_date"])

    return data


def compute_publication_production(data: pd.DataFrame):
    """
    Compute monthly and yearly counts of publications for each 'pi_id' and 'seed'.

    Args:
        data (pd.DataFrame): The input DataFrame containing publication data.

    Returns:
        tuple: A tuple containing two DataFrames - monthly_publications and
            yearly_publications.
            - monthly_publications
            - yearly_publications
    """

    # Set 'publication_date' as the index
    data.set_index("publication_date", inplace=True)

    # Group by 'pi_id' and 'seed', then resample to get monthly and yearly counts of publications
    grouped = data.groupby(["pi_id", "seed"])
    monthly_publications = grouped.resample("M").size().reset_index(name="count")
    yearly_publications = grouped.resample("Y").size().reset_index(name="count")

    return monthly_publications, yearly_publications


def get_event_study_outputs(data: pd.DataFrame, level0: pd.DataFrame, output_type: str):
    """
    Compute event study outputs based on the provided data and level0 data.

    Args:
        data (pd.DataFrame): The input data.
        level0 (pd.DataFrame): The level0 data.

    Returns:
        pd.DataFrame: The computed event study outputs.

    """
    # convert level0 publication_date and parent_publication_date to datetime
    level0["publication_date"] = pd.to_datetime(level0["publication_date"])
    level0["parent_publication_date"] = pd.to_datetime(
        level0["parent_publication_date"]
    )

    logger.info("Merging data with level0 data")
    merged_data = pd.merge(
        data,
        level0[["id", "parent_id", "parent_publication_date"]],
        on="id",
        how="inner",
    )

    logger.info("Sorting data by 'parent_id' and 'parent_publication_date'")
    merged_data.sort_values(
        by=["parent_id", "parent_publication_date"],
        ascending=[True, True],
        inplace=True,
    )

    logger.info(
        "Computing time difference between 'publication_date' and 'parent_publication_date'"
    )
    parent_id_in_list = merged_data["parent_id"].isin(
        ["W3177828909", "W3211795435", "W3202105508"]
    )

    logger.info("Reordering data based on 'parent_id'")
    merged_data = pd.concat(
        [merged_data.loc[parent_id_in_list], merged_data.loc[~parent_id_in_list]]
    )

    logger.info("Dropping duplicates based on 'pi_id'")
    merged_data.drop_duplicates(subset="pi_id", keep="first", inplace=True)

    logger.info("Mapping 'parent_publication_date' to 'pi_id'")
    pi_id_date_dict = dict(
        zip(merged_data["pi_id"], merged_data["parent_publication_date"])
    )
    data["parent_publication_date"] = data["pi_id"].map(pi_id_date_dict)

    logger.info(
        "Computing time difference between 'publication_date' and 'parent_publication_date'"
    )
    data["time"] = (
        data["publication_date"] - data["parent_publication_date"]
    ) / pd.Timedelta(days=90)

    # Convert to integer
    final_data = data.copy()
    final_data.dropna(subset=["time"], inplace=True)
    final_data["time"] = final_data["time"].astype(int)
    if output_type == "publications":
        final_data = (
            final_data.groupby(["pi_id", "time", "seed"])
            .size()
            .reset_index(name="count")
        )
    elif output_type == "citations":
        final_data = (
            final_data.groupby(["pi_id", "time", "seed"])["cited_by_count"]
            .sum()
            .reset_index(name="count")
        )

    return data, final_data


def get_event_study_strength(data, sc_data_af, sc_data_ct, output_type):
    """
    Compute event study strength based on the provided data and strong links data.

    Args:
        data (pd.DataFrame): The input data.
        sc_data_af (pd.DataFrame): The strong links data for AlphaFold.
        sc_data_ct (pd.DataFrame): The strong links data for counterfactual papers.
        output_type (str): The output type.

    Returns:
        pd.DataFrame: The computed event study strength.

    """
    sc_data = pd.concat([sc_data_af, sc_data_ct])

    # drop none in strength
    sc_data = sc_data.dropna(subset=["strength"])
    # drop duplicates
    sc_data = sc_data.drop_duplicates(subset=["id"], keep="first")

    data["strong"] = data["id"].isin(sc_data["id"])

    # for PI with at least a strong link, set all papers to strong
    pi_strong = data[data["strong"]]["pi_id"].unique()
    data.loc[data["pi_id"].isin(pi_strong), "strong"] = True

    if output_type == "publications":
        final_data = (
            data.groupby(["pi_id", "time", "seed", "strong"])
            .size()
            .reset_index(name="count")
        )
    elif output_type == "citations":
        final_data = (
            data.groupby(["pi_id", "time", "seed", "strong"])["cited_by_count"]
            .sum()
            .reset_index(name="count")
        )

    return final_data


def get_event_study_pdb_submissions(
    data: pd.DataFrame, pdb_data: pd.DataFrame, output_type: str
):
    """
    Compute event study outputs based on the provided data and pdb_data.

    Args:
        data (pd.DataFrame): The input data.
        pdb_data (pd.DataFrame): The pdb data.
        output_type (str): The output type.

    Returns:
        pd.DataFrame: The computed event study outputs.
    """

    logger.info("Merging data with pdb_data data")

    data = pd.merge(
        data,
        pdb_data[["id", "resolution", "R_free"]],
        on="id",
        how="inner",
    )

    if output_type == "publications":
        final_data = (
            data.groupby(["pi_id", "time", "seed"]).size().reset_index(name="count")
        )
    elif output_type == "citations":
        final_data = (
            data.groupby(["pi_id", "time", "seed"])["cited_by_count"]
            .sum()
            .reset_index(name="count")
        )
    return final_data
