"""
This is a boilerplate pipeline 'data_processing_chains'
generated using Kedro 0.19.1
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import numpy as np
from Bio import Entrez
from .utils import (
    sort_and_drop,
    build_chain_dict,
    flatten_dict,
    ensure_levels,
    breaks_chain,
    transform_long,
    transform_long_pairs,
)

pd.options.mode.copy_on_write = True

Entrez.email = "david.ampudia@nesta.org.uk"
logger = logging.getLogger(__name__)


def _process_seed_id(
    seed_id: str, depth_data: pd.DataFrame, identifier: str, num_levels: int
) -> pd.DataFrame:
    """
    Worker function to process a single seed_id.

    Args:
        seed_id (str): The seed identifier to process.
        depth_data (pd.DataFrame): The input DataFrame containing citation link data.
        identifier (str): The identifier used to filter the citation links.
        num_levels (int): The number of levels in the citation chain.

    Returns:
        pd.DataFrame: The processed DataFrame for this seed_id.
    """
    logger.info("Building citation chain for seed %s", seed_id)
    data_dict = build_chain_dict(depth_data, identifier, seed_id, 0, num_levels)
    flat_rows = flatten_dict(data_dict)
    seed_df = pd.DataFrame(flat_rows)
    seed_df["seed_id"] = seed_id
    seed_df = ensure_levels(seed_df, num_levels)
    seed_df = seed_df[
        [
            "seed_id",
        ]
        + [f"level_{i}" for i in range(num_levels)]
    ]
    return seed_df


def filter_relevant_citation_links(
    depth_data: pd.DataFrame,
    identifier: str,
    num_levels: int,
    max_workers: int = 8,
    **kwargs,
) -> pd.DataFrame:
    """
    Filters and processes citation links based on relevance.

    Args:
        depth_data (pd.DataFrame): The input DataFrame containing citation
             link data.
        identifier (str): The identifier used to filter the citation links.
        num_levels (int): The number of levels in the citation chain.
        max_workers (int): Maximum number of threads to use for parallel processing.
            Defaults to 8.
        **kwargs: Additional keyword arguments for custom sorting.

    Returns:
        pd.DataFrame: The processed DataFrame containing relevant citation chains.
    """

    logger.info("Exploding citation links")
    depth_data = depth_data.explode("strength")

    # logger.info("Remove citation links with no strength")
    # depth_data.dropna(subset=["strength"], inplace=True)

    for col in ["intent"]:
        logger.info("Extracting %s from citation links", col)
        depth_data[col] = depth_data["strength"].apply(
            lambda x: x[col] if x else None  # pylint: disable=cell-var-from-loop
        )

    # create strong variable
    depth_data["intent"] = depth_data["intent"].apply(
        lambda x: (
            "strong"
            if x in ["result", "methodology"]
            else (
                "weak" if x in ["background"] else "unknown" if x in [""] else "no_data"
            )
        )
    )

    # define & implement the custom sorting order
    depth_data = sort_and_drop(depth_data, **kwargs)

    if identifier == "pmid":
        # assign 99999999 PMID to Multimer paper, with OA id W3202105508
        depth_data.loc[depth_data["parent_id"] == "W3202105508", "parent_pmid"] = (
            "99999999"
        )

    logger.info("Filter nulls for relevant identifier %s", identifier)
    depth_data = depth_data[depth_data[identifier].notnull()]

    logger.info("Drop duplicates of parent_id, identifier, and intent")
    depth_data.drop_duplicates(subset=["parent_id", identifier, "intent"], inplace=True)

    seeds = depth_data[depth_data["level"] == 0][f"parent_{identifier}"].unique()

    if len(seeds) == 0:
        # move level up by 1
        depth_data["level"] = depth_data["level"] - 1

    # get the final list of seeds after potential level adjustment
    seed_list = list(
        depth_data[depth_data["level"] == 0][f"parent_{identifier}"].unique()
    )

    logger.info("Processing %d seeds using multithreading", len(seed_list))
    output = []
    max_workers = min(len(seed_list), max_workers)  # limit to max_workers threads
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # submit all tasks
        future_to_seed = {
            executor.submit(
                _process_seed_id, seed_id, depth_data, identifier, num_levels
            ): seed_id
            for seed_id in seed_list
        }

        # collect results as they complete
        for future in as_completed(future_to_seed):
            seed_id = future_to_seed[future]
            try:
                result = future.result()
                output.append(result)
                logger.info("Completed processing seed %s", seed_id)
            except Exception as exc:
                logger.error("Seed %s generated an exception: %s", seed_id, exc)
                raise

    chains_df = pd.concat(output)

    # for each level, merge chains with paper intent data
    for i in range(num_levels):
        if i == 0:
            chains_df = chains_df.merge(
                depth_data[[f"parent_{identifier}", identifier, "intent"]],
                left_on=["seed_id", "level_0"],
                right_on=[f"parent_{identifier}", identifier],
                how="left",
            )
        else:
            chains_df = chains_df.merge(
                depth_data[[f"parent_{identifier}", identifier, "intent"]],
                left_on=[f"level_{i-1}", f"level_{i}"],
                right_on=[f"parent_{identifier}", identifier],
                how="left",
            )
        # rename the intent column
        chains_df.rename(columns={"intent": f"intent_{i}"}, inplace=True)

        # drop duplicate columns
        chains_df.drop(columns=[f"parent_{identifier}", identifier], inplace=True)

    # substitute empty string in intent for None
    for i in range(num_levels):
        chains_df[f"intent_{i}"] = chains_df[f"intent_{i}"].apply(
            lambda x: x if x != "" else None
        )

    # in the intent columns, replace NaN with 'N/A'
    for i in range(num_levels):
        chains_df[f"intent_{i}"] = chains_df[f"intent_{i}"].apply(
            lambda x: "N/A" if x is np.nan else x
        )
    complete_chains_df = chains_df[
        ~chains_df.apply(lambda row: breaks_chain(row, num_levels), axis=1)
    ]

    # drop rows that have all four intent as nan or N/A
    columns_to_check = (
        ["intent_0", "intent_1", "intent_2"]
        if num_levels > 2
        else ["intent_0", "intent_1"]
    )

    complete_chains_df = complete_chains_df[
        ~complete_chains_df[columns_to_check].apply(
            lambda row: all(map(lambda x: x == "N/A" or pd.isna(x), row)), axis=1
        )
    ]

    return complete_chains_df


def get_papers_with_strong_chain(
    chains: pd.DataFrame,
    depth_data: pd.DataFrame,
    identifier: str,
    num_levels: int,
) -> pd.DataFrame:
    """
    Retrieves papers with the full chain of data for a given identifier.

    Args:
        chains (pd.DataFrame): DataFrame containing the chains data.
        depth_data (pd.DataFrame): DataFrame containing the AlphaFold data.
        identifier (str): Identifier used to match the chains and AlphaFold data.
        num_levels (int): Number of levels in the chains DataFrame.

    Returns:
        pd.DataFrame: DataFrame containing the merged data of chains and AlphaFold data.
    """
    logger.info("Transforming long chains for identifier %s", identifier)
    data_ids = transform_long(chains, ["strong"], num_levels)

    logger.info("Merging chains with AlphaFold data")
    depth_data = depth_data.merge(
        data_ids,
        how="inner",
        left_on=[identifier, "level"],
        right_on=["paper_id", "level"],
    )
    return depth_data


def get_chain_label_papers(
    chains: pd.DataFrame,
    depth_data: pd.DataFrame,
    identifier: str,
    num_levels: int,
) -> pd.DataFrame:
    """
    Identifies and labels chains in the depth_data DataFrame based on the information
    in the chains DataFrame.

    Args:
        chains (pd.DataFrame): DataFrame containing information about chains.
        depth_data (pd.DataFrame): DataFrame containing depth data.
        identifier (str): Identifier column name in the depth_data DataFrame.
        num_levels (int): Number of levels in the chains DataFrame.

    Returns:
        pd.DataFrame: DataFrame with the labeled chains.

    Raises:
        None
    """
    logger.info("Identifying strong-only chains")
    data_ids_strong = transform_long_pairs(chains, ["strong"], num_levels)

    logger.info("Identifying weak-only chains")
    data_ids_weak = transform_long_pairs(chains, ["weak"], num_levels)

    logger.info("Identifying unknown only chains")
    data_ids_unknown = transform_long_pairs(chains, ["unknown"], num_levels)

    logger.info("Identifying no_data chains only")
    data_ids_no_data = transform_long_pairs(chains, ["no_data"], num_levels)

    logger.info("Identifying partial only-strong chains")
    data_ids_partial_strong = transform_long_pairs(
        chains, ["strong", "no_data", "unknown"], num_levels
    )

    # antijoin
    data_ids_partial_strong = data_ids_partial_strong[
        ~(
            (
                (data_ids_partial_strong["paper_id"].isin(data_ids_strong["paper_id"]))
                & (
                    data_ids_partial_strong["parent_paper_id"].isin(
                        data_ids_strong["parent_paper_id"]
                    )
                )
            )
            | (
                (data_ids_partial_strong["paper_id"].isin(data_ids_no_data["paper_id"]))
                & (
                    data_ids_partial_strong["parent_paper_id"].isin(
                        data_ids_no_data["parent_paper_id"]
                    )
                )
            )
            | (
                (data_ids_partial_strong["paper_id"].isin(data_ids_unknown["paper_id"]))
                & (
                    data_ids_partial_strong["parent_paper_id"].isin(
                        data_ids_unknown["parent_paper_id"]
                    )
                )
            )
        )
    ]

    logger.info("Identify partial weak chains")
    data_ids_partial_weak = transform_long_pairs(
        chains, ["weak", "no_data", "unknown"], num_levels
    )

    # antijoin
    data_ids_partial_weak = data_ids_partial_weak[
        ~(
            (
                (data_ids_partial_weak["paper_id"].isin(data_ids_weak["paper_id"]))
                & (
                    data_ids_partial_weak["parent_paper_id"].isin(
                        data_ids_weak["parent_paper_id"]
                    )
                )
            )
            | (
                (data_ids_partial_weak["paper_id"].isin(data_ids_no_data["paper_id"]))
                & (
                    data_ids_partial_weak["parent_paper_id"].isin(
                        data_ids_no_data["parent_paper_id"]
                    )
                )
            )
            | (
                (data_ids_partial_weak["paper_id"].isin(data_ids_unknown["paper_id"]))
                & (
                    data_ids_partial_weak["parent_paper_id"].isin(
                        data_ids_unknown["parent_paper_id"]
                    )
                )
            )
        )
    ]

    logger.info("Identifying partial-strong chains")
    data_ids_mixed = transform_long_pairs(chains, ["strong", "weak"], num_levels)

    # antijoin
    data_ids_mixed = data_ids_mixed[
        ~(
            (
                data_ids_mixed["paper_id"].isin(data_ids_strong["paper_id"])
                & data_ids_mixed["parent_paper_id"].isin(
                    data_ids_strong["parent_paper_id"]
                )
            )
            | (
                data_ids_mixed["paper_id"].isin(data_ids_weak["paper_id"])
                & data_ids_mixed["parent_paper_id"].isin(
                    data_ids_weak["parent_paper_id"]
                )
            )
        )
    ]

    # add a column in each
    data_ids_strong["chain_label"] = "strong"
    data_ids_weak["chain_label"] = "weak"
    data_ids_unknown["chain_label"] = "unknown"
    data_ids_no_data["chain_label"] = "no_data"
    data_ids_partial_strong["chain_label"] = "partial_strong"
    data_ids_partial_weak["chain_label"] = "partial_weak"
    data_ids_mixed["chain_label"] = "mixed"

    # concatenate all
    chain_labels = pd.concat(
        [
            data_ids_strong,
            data_ids_weak,
            data_ids_unknown,
            data_ids_no_data,
            data_ids_partial_strong,
            data_ids_partial_weak,
            data_ids_mixed,
        ]
    )

    # map
    chain_labels_map = chain_labels.set_index(["paper_id", "parent_paper_id"])[
        "chain_label"
    ].to_dict()

    # map chain_label to depth_data identifier, parent_{identifier}
    depth_data["chain_label"] = depth_data.apply(
        lambda x: chain_labels_map.get(
            (x[identifier], x[f"parent_{identifier}"]), None
        ),
        axis=1,
    )

    # fillna chain_label with "other"
    depth_data["chain_label"].fillna("other", inplace=True)

    return depth_data
