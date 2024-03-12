"""
This is a boilerplate pipeline 'data_analysis_exploratory'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict, Callable
import pandas as pd
from kedro.io import AbstractDataset

logger = logging.getLogger(__name__)


def combine_depth_strength_level_0(
    depth_data: pd.DataFrame,
    strength_data: Dict[str, Callable[[], pd.DataFrame]],
) -> pd.DataFrame:
    """
    Combines depth data with strength data at level 0.

    Args:
        depth_data (pd.DataFrame): The depth data to be combined.
        strength_data (Dict[str, Callable[[], pd.DataFrame]]): The strength data to be combined.

    Returns:
        pd.DataFrame: The merged dataframe containing the combined data.

    """
    logger.info("Combining depth and strength data")
    depth_data["parent_doi"] = "10.1038/s41586-021-03819-2"
    strength_data_grouped = (
        strength_data.groupby(["parent_doi", "child_doi"])[
            ["influential", "intent", "context"]
        ]
        .apply(lambda x: x.to_dict(orient="records"))
        .reset_index()
        .rename(columns={0: "strength"})
    )

    # drop parent_doi
    strength_data_grouped = strength_data_grouped.drop(["parent_doi"], axis=1)

    merged_df = pd.merge(
        depth_data,
        strength_data_grouped,
        left_on="doi",
        right_on="child_doi",
        how="left",
    )
    merged_df = merged_df.drop(["child_doi"], axis=1)

    logger.info("Merging depth and strength data complete")
    return merged_df[
        [
            "parent_id",
            "id",
            "parent_doi",
            "doi",
            "level",
            "publication_date",
            "mesh_terms",
            "strength",
            "cited_by_count"
        ]
    ]


def combine_depth_strength_other_levels(
    previous_depth_data: pd.DataFrame,
    depth_data: pd.DataFrame,
    strength_data: Dict[str, Callable[[], pd.DataFrame]],
) -> pd.DataFrame:
    """
    Combines depth data, strength data, and other level data into a single DataFrame.

    Args:
        previous_depth_data (pd.DataFrame): DataFrame containing previous depth data.
        depth_data (pd.DataFrame): DataFrame containing depth data.
        strength_data (Dict[str, Callable[[], pd.DataFrame]]): Dictionary of strength data.

    Returns:
        pd.DataFrame: Merged DataFrame containing combined data.

    """
    logger.info("Combining depth and strength data")
    previous_depth_data_unique = previous_depth_data.drop_duplicates(subset=["id"])
    depth_data["parent_doi"] = depth_data["parent_id"].map(
        previous_depth_data_unique.set_index("id")["doi"]
    )

    logger.info("Loading strength data. Total number of works: %s", len(strength_data))
    strength_data = pd.concat([v() for v in strength_data.values()])

    strength_data_grouped = (
        strength_data.groupby(["parent_doi", "child_doi"])[
            ["influential", "intent", "context"]
        ]
        .apply(lambda x: x.to_dict(orient="records"))
        .reset_index()
        .rename(columns={0: "strength"})
    )
    merged_df = pd.merge(
        depth_data,
        strength_data_grouped,
        left_on=["parent_doi", "doi"],
        right_on=["parent_doi", "child_doi"],
        how="left",
    )
    merged_df = merged_df.drop(["child_doi"], axis=1)
    logger.info("Merging depth and strength data complete")
    return merged_df[
        [
            "parent_id",
            "id",
            "parent_doi",
            "doi",
            "level",
            "publication_date",
            "mesh_terms",
            "strength",
            "cited_by_count"
        ]
    ]


def process_subfield_data(data: Dict[str, AbstractDataset]) -> pd.DataFrame:
    """
    Process data by level without enriching the data with additional mesh terms
        and return a DataFrame with selected columns.

    Args:
        data (Dict[str, AbstractDataset]): A dictionary containing the input data.
        level (int): The level to process the data for.

    Returns:
        pd.DataFrame: A DataFrame containing the processed data with selected columns.
    """

    logger.info("Processing subfield data. Total number of laoders: %s", len(data))
    # Initialize an empty DataFrame for current level
    subfield_df = pd.DataFrame()

    for i, loader in enumerate(data.values()):
        logger.info("Processing subfield data. Loader: %s / %s", i + 1, len(data))
        raw_json_data = loader()
        json_data = [
            {
                k: v
                for k, v in item.items()
                if k in ["id", "doi", "publication_date", "mesh_terms", "cited_by_count"]
            }
            for item in raw_json_data
        ]

        # transform to dataframe and add parent_id
        df = pd.DataFrame(json_data)

        # check mesh_terms exists else skip
        if "mesh_terms" not in df.columns:
            logger.warning("Skipping subfield data processing. No mesh terms found.")
            continue

        # mesh tuples
        df["mesh_terms"] = df["mesh_terms"].apply(
            lambda x: (
                [(c["descriptor_ui"], c["descriptor_name"]) for c in x] if x else None
            )
        )

        # change doi to remove the url
        df["doi"] = df["doi"].str.replace("https://doi.org/", "")

        # change id to remove openalex url
        df["id"] = df["id"].str.replace("https://openalex.org/", "")
        

        # Append current batch to level DataFrame
        subfield_df = pd.concat([subfield_df, df])

    return subfield_df[["id", "doi", "publication_date", "mesh_terms", "cited_by_count"]]
