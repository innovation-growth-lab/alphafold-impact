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
    process_intent_references: Process the reference outputs with simplified intents field and generate a list of rows
        containing relevant information.
"""

import logging
import asyncio
import pandas as pd
from typing import Dict
from kedro.io import AbstractDataset
from .utils import (
    get_intent_level_0_async,
    get_intent_level_n_async,
)

logger = logging.getLogger(__name__)


def get_citation_intent_from_oa_dataset(
    oa_dataset: pd.DataFrame,
    **kwargs,
) -> pd.DataFrame:
    """
    Retrieves citation intent data from an Open Access dataset.

    Args:
        oa_dataset (pd.DataFrame): The Open Access dataset containing citation
            information.
        **kwargs: Additional keyword arguments.

    Returns:
        pd.DataFrame: The processed dataset with citation intent information.
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

    # create a dictionary to map parent_id to DOI and PMID
    parent_info = {
        "W3177828909": {"doi": "10.1038/s41586-021-03819-2", "pmid": "34265844"},
        "W3211795435": {"doi": "10.1093/nar/gkab1061", "pmid": "34791371"},
        "W3202105508": {"doi": "10.1101/2021.10.04.463034", "pmid": None},
    }

    # replace the values in the 'parent_doi' and 'parent_pmid' columns when level is 0
    oa_dataset.loc[oa_dataset["level"] == 0, "parent_doi"] = oa_dataset.loc[
        oa_dataset["level"] == 0, "parent_id"
    ].map(lambda x: parent_info[str(x)]["doi"])
    oa_dataset.loc[oa_dataset["level"] == 0, "parent_pmid"] = oa_dataset.loc[
        oa_dataset["level"] == 0, "parent_id"
    ].map(lambda x: parent_info[str(x)]["pmid"])

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

    # Read all level data
    level_dfs = {level: loader() for level, loader in levels.items()}

    # Concatenate all levels
    processed_df = pd.concat([df for df in level_dfs.values()])

    # Log empty DOIs with non-empty PMIDs
    logger.info(
        "Number of rows with empty doi and non-empty pmid: %d",
        processed_df[
            (processed_df["doi"] == None) & (processed_df["pmid"] != None)
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


def get_baseline_seed_intent(oa_dataset: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """
    Retrieves the seed intent data from the given baseline OA dataset.

    Args:
        oa_dataset (pd.DataFrame): The input OA dataset.
        **kwargs: Additional keyword arguments.

    Returns:
        pd.DataFrame: The processed seed intent data.
    """
    return None

    # # remove duplicate rows based on "id"
    # oa_dataset = oa_dataset.drop_duplicates(subset="id")

    # inputs = oa_dataset.apply(
    #     lambda x: (x["id"], "", x["doi"], "", x["pmid"]),
    #     axis=1,
    # ).tolist()

    # # use joblib to parallelize the function calls
    # level_outputs = Parallel(n_jobs=4)(
    #     delayed(iterate_citation_detail_points)(*input, direction="citations", **kwargs)
    #     for input in inputs
    # )

    # # change oa_dataset "doi" to "parent_doi"
    # oa_dataset.rename(columns={"doi": "parent_doi"}, inplace=True)
    # level_dict = dict(list(zip(oa_dataset["parent_doi"], level_outputs)))

    # processed_level_citations = process_citations(level_dict)

    # return pd.DataFrame(
    #     processed_level_citations,
    #     columns=[
    #         "parent_doi",
    #         "pmid",
    #         "doi",
    #         "influential",
    #         "intent",
    #         "context",
    #     ],
    # )


def get_baseline_citation_intent_from_oa_dataset(
    oa_dataset: pd.DataFrame,
    **kwargs,
) -> pd.DataFrame:
    """
    Retrieves citation intent data from an Open Access dataset.

    Args:
        oa_dataset (pd.DataFrame): The Open Access dataset containing citation
            information.
        **kwargs: Additional keyword arguments.

    Returns:
        pd.DataFrame: The processed dataset with citation intent information.

    """
    return None
    # # if level is str, convert to int with "seed" as -1
    # oa_dataset["level"] = oa_dataset["level"].apply(
    #     lambda x: -1 if x == "seed" else int(x)
    # )

    # # Create a mapping from id to doi and pmid for each level
    # oa_dataset["parent_level"] = oa_dataset["level"] - 1

    # oa_dataset = pd.merge(
    #     oa_dataset,
    #     oa_dataset,
    #     left_on=["parent_id", "parent_level"],
    #     right_on=["id", "level"],
    #     how="left",
    #     suffixes=("", "_parent"),
    # )
    # oa_dataset.rename(
    #     columns={"doi_parent": "parent_doi", "pmid_parent": "parent_pmid"}, inplace=True
    # )

    # # drop the other _parent columns
    # oa_dataset.drop(
    #     columns=[col for col in oa_dataset.columns if "_parent" in col], inplace=True
    # )

    # # get the citation intent for each level
    # processed_df = get_intent(oa_dataset, **kwargs)

    # # check how often doi is empty and pmid is not
    # logger.info(
    #     "Number of rows with empty doi and non-empty pmid: %d",
    #     processed_df[(processed_df["doi"] == "") & (processed_df["pmid"] != "")].shape[
    #         0
    #     ],
    # )

    # processed_grouped_data_doi = (
    #     processed_df.groupby(["parent_doi", "doi"])[
    #         ["influential", "intent", "context"]
    #     ]
    #     .apply(lambda x: x.to_dict(orient="records"))
    #     .reset_index()
    #     .rename(columns={0: "strength"})
    # )

    # processed_grouped_data_pmid = (
    #     processed_df.groupby(["parent_doi", "pmid"])[
    #         ["influential", "intent", "context"]
    #     ]
    #     .apply(lambda x: x.to_dict(orient="records"))
    #     .reset_index()
    #     .rename(columns={0: "strength"})
    # )

    # processed_data = pd.merge(
    #     oa_dataset, processed_grouped_data_doi, on=["parent_doi", "doi"], how="left"
    # )

    # processed_data = pd.merge(
    #     processed_data,
    #     processed_grouped_data_pmid,
    #     left_on=["parent_doi", "parent_pmid"],
    #     right_on=["parent_doi", "pmid"],
    #     how="left",
    # )

    # # combine the two strength columns
    # processed_data["strength"] = processed_data["strength_x"].combine_first(
    #     processed_data["strength_y"]
    # )

    # # rename the pmid columns
    # processed_data.rename(columns={"pmid_x": "pmid"}, inplace=True)
    # processed_data.drop(columns=["strength_x", "strength_y", "pmid_y"], inplace=True)

    # return processed_data[
    #     [
    #         "parent_id",
    #         "parent_doi",
    #         "parent_pmid",
    #         "id",
    #         "doi",
    #         "pmid",
    #         "level",
    #         "publication_date",
    #         "mesh_terms",
    #         "cited_by_count",
    #         "authorships",
    #         "parent_level",
    #         "strength",
    #         "topics",
    #         "concepts",
    #     ]
    # ]
