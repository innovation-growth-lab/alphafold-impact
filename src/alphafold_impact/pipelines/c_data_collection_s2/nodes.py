"""Pipeline nodes for the data_collection_s2 pipeline.

Functions:
    fetch_citation_details: Fetches citation details for a given work ID.
    process_citations: Process citation outputs and extract relevant
        information.
    process_references: Process the reference outputs and generate a
        list of rows containing relevant information.
    iterate_citation_detail_points: Iterates over citation detail
        points and fetches citation details based on the provided parameters.
    get_intent_level_0: Retrieves the intent level 0 data from the given
        OA dataset.
    get_intent_level: Retrieves the intent level data from the given
        OA dataset based on the specified level.
    get_citation_intent_from_oa_dataset: Retrieves citation intent data
        from an Open Access dataset.
    process_intent_references: Process the reference outputs with simplified intents
        field and generate a list of rows containing relevant information.
    process_citation_levels: Process and combine citation data from all levels.
    get_baseline_seed_intent: Retrieves the seed intent data from the given baseline OA dataset.
    get_baseline_intent: Retrieves citation intent data for other baseline papers.
"""

import logging
import asyncio
from typing import Dict, Generator
import pandas as pd
from kedro.io import AbstractDataset
from .utils import (
    create_parent_child_df,
    get_intent_level_0_async,
    get_intent_level_n_async,
)

logger = logging.getLogger(__name__)


def get_citation_intent_from_oa_dataset(
    oa_dataset: pd.DataFrame,
    **kwargs,
) -> Generator[Dict[str, AbstractDataset], None, None]:
    """
    Retrieves citation intent data from an Open Access dataset.

    Args:
        oa_dataset (pd.DataFrame): The Open Access dataset containing citation
            information.
        **kwargs: Additional keyword arguments.

    Returns:
        Generator[Dict[str, AbstractDataset], None, None]: Generator containing the
            level data for partitioned storage.
    """
    oa_dataset = create_parent_child_df(oa_dataset)

    # Create task definitions for all levels
    task_definitions = [
        lambda: get_intent_level_0_async(oa_dataset, **kwargs),
        lambda: get_intent_level_n_async(oa_dataset, 1, **kwargs),
        lambda: get_intent_level_n_async(oa_dataset, 2, **kwargs),
    ]

    # Process tasks as they complete
    for i, task_def in enumerate(task_definitions):
        level_name = f"level_{i}"
        level_df = asyncio.run(task_def())
        level_df = level_df.replace("", None)
        yield {level_name: level_df}


def process_citation_levels(
    oa_dataset: pd.DataFrame, levels: Dict[str, AbstractDataset]
) -> pd.DataFrame:
    """
    Process and combine citation data from all levels.

    Args:
        oa_dataset: Original dataset with paper information
        level_paths: Dictionary containing paths to saved level data files

    Returns:
        pd.DataFrame: Processed and combined citation data
    """

    oa_dataset = create_parent_child_df(oa_dataset)

    # Read all level data
    level_dfs = {level: loader() for level, loader in levels.items()}

    # Concatenate all levels
    processed_df = pd.concat([df for df in level_dfs.values()])

    # Log empty DOIs with non-empty PMIDs
    logger.info(
        "Number of rows with empty doi and non-empty pmid: %d",
        processed_df[
            (processed_df["doi"] == None)  # pylint: disable=C0121
            & (processed_df["pmid"] != None)  # pylint: disable=C0121
        ].shape[0],
    )

    # Process DOI-based strengths
    processed_grouped_data_doi = (
        processed_df.groupby(["parent_doi", "doi"])[["influential", "intent"]]
        .apply(lambda x: x.to_dict(orient="records"))
        .reset_index()
        .rename(columns={0: "strength"})
    )

    # Process PMID-based strengths
    processed_grouped_data_pmid = (
        processed_df.groupby(["parent_doi", "pmid"])[["influential", "intent"]]
        .apply(lambda x: x.to_dict(orient="records"))
        .reset_index()
        .rename(columns={0: "strength"})
    )

    # Merge with original dataset
    processed_data = pd.merge(
        oa_dataset, processed_grouped_data_doi, on=["parent_doi", "doi"], how="left"
    )

    processed_data = pd.merge(
        processed_data,
        processed_grouped_data_pmid,
        left_on=["parent_doi", "parent_pmid"],
        right_on=["parent_doi", "pmid"],
        how="left",
    )

    # Combine strength columns
    processed_data["strength"] = processed_data["strength_x"].combine_first(
        processed_data["strength_y"]
    )

    # Clean up columns
    processed_data.rename(columns={"pmid_x": "pmid"}, inplace=True)
    processed_data.drop(columns=["strength_x", "strength_y", "pmid_y"], inplace=True)

    return processed_data[
        [
            "parent_id",
            "parent_doi",
            "parent_pmid",
            "id",
            "doi",
            "pmid",
            "level",
            "publication_date",
            "mesh_terms",
            "cited_by_count",
            "authorships",
            "parent_level",
            "strength",
            "topics",
            "concepts",
        ]
    ]


def get_baseline_seed_intent(
    oa_dataset: pd.DataFrame, **kwargs
) -> Generator[AbstractDataset, None, None]:
    """
    Retrieves the seed intent data from the given baseline OA dataset.
    Uses async processing and yields level results for partitioned storage.

    Args:
        oa_dataset (pd.DataFrame): The input OA dataset.
        **kwargs: Additional keyword arguments.

    Yields:
        Generator[AbstractDataset, None, None]: Generator containing the level data.
    """
    # remove duplicate rows based on "id"
    oa_dataset = oa_dataset.drop_duplicates(subset="id")

    # remove the triad of AF DOIs
    oa_dataset = oa_dataset[
        ~oa_dataset["doi"].isin(
            [
                "10.1038/s41586-021-03819-2",
                "10.1093/nar/gkab1061",
                "10.1101/2021.10.04.463034",
            ]
        )
    ]

    # rename id to parent_id, pmid to parent_pmid, doi to parent_doi
    oa_dataset.rename(
        columns={"id": "parent_id", "pmid": "parent_pmid", "doi": "parent_doi"},
        inplace=True,
    )

    level_df = asyncio.run(get_intent_level_n_async(oa_dataset, None, **kwargs))
    level_df = level_df.replace("", None)
    return level_df


def get_baseline_intent(
    oa_dataset: pd.DataFrame, **kwargs
) -> Generator[AbstractDataset, None, None]:
    """
    Retrieves citation intent data for other baseline papers.
    Uses async processing and yields level results for partitioned storage.

    Args:
        oa_dataset (pd.DataFrame): The input OA dataset.
        **kwargs: Additional keyword arguments.

    Yields:
        Generator[AbstractDataset, None, None]: Generator containing the level data.
    """
    # if level is str, convert to int with "seed" as -1
    oa_dataset["level"] = oa_dataset["level"].apply(
        lambda x: -1 if x == "seed" else int(x)
    )

    # Create a mapping from id to doi and pmid for each level
    oa_dataset["parent_level"] = oa_dataset["level"] - 1

    oa_dataset = pd.merge(
        oa_dataset,
        oa_dataset,
        left_on=["parent_id", "parent_level"],
        right_on=["id", "level"],
        how="left",
        suffixes=("", "_parent"),
    )
    oa_dataset.rename(
        columns={"doi_parent": "parent_doi", "pmid_parent": "parent_pmid"}, inplace=True
    )

    # drop the other _parent columns
    oa_dataset.drop(
        columns=[col for col in oa_dataset.columns if "_parent" in col], inplace=True
    )

    # Create task definitions for all levels
    task_definitions = {
        "level_1": lambda: get_intent_level_n_async(oa_dataset, 1, **kwargs),
        "level_2": lambda: get_intent_level_n_async(oa_dataset, 2, **kwargs),
    }

    # Process tasks as they complete
    for level_name, task_def in task_definitions.items():
        level_df = asyncio.run(task_def())
        level_df = level_df.replace("", None)
        yield {level_name: level_df}
