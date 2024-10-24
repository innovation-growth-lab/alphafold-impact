"""
This module provides utility functions for data processing chains in the 
AlphaFold Impact project.
"""

import logging
from typing import Dict
import pandas as pd
import numpy as np

pd.options.mode.copy_on_write = True
logger = logging.getLogger(__name__)


def sort_and_drop(
    data: pd.DataFrame,
    unique: bool = False,
):
    """
    Sorts the data based on a custom sorting order and drops duplicate citation links if specified.

    Args:
        data (pd.DataFrame): The input DataFrame to be sorted and processed.
        unique (bool, optional): If True, drops duplicate citation links. Defaults to False.

    Returns:
        pd.DataFrame: The sorted and processed DataFrame.

    """
    # define & implement the custom sorting order
    # sort_order = {"methodology": 0, "result": 1, "background": 2, "": 3, None: 4}
    sort_order = {"strong": 0, "weak": 1, "unknown": 2, "N/A": 3}
    data["sort_order"] = data["intent"].map(sort_order)

    data = (
        data.sort_values("sort_order")
        .groupby(["parent_id", "id"])
        .first()
        .reset_index()
    )

    # drop W31
    # drop rows with id W3177828909, W3211795435, W3202105508
    data = data[
        ~data["id"].isin(["W3177828909", "W3211795435", "W3202105508", "W3202105508"])
    ]

    # drop rows with parent_id W3177828909, W3211795435, W3202105508 except in level 0
    data = data[
        ~(
            (data["parent_id"].isin(["W3177828909", "W3211795435", "W3202105508"]))
            & (data["level"] != 0)
        )
    ]

    if unique:
        logger.info("Dropping duplicate citation links")
        data.sort_values("sort_order").drop_duplicates(subset="id", inplace=True)

    # drop the sort_order column, reindex
    data.drop(columns="sort_order", inplace=True)
    data.reset_index(drop=True, inplace=True)

    return data


def build_chain_dict(df, identifier, parent, level, num_levels):
    """
    Recursively builds a dictionary representing a chain of data processing nodes.

    Args:
        df (pandas.DataFrame): The DataFrame containing the data.
        identifier (str): The column name representing the unique identifier for each node.
        parent (str): The parent identifier for the current level.
        level (int): The current level of the chain.

    Returns:
        dict: A dictionary representing the chain of data processing nodes.
    """
    # base case: if level is greater than 2, return an empty dictionary
    if level >= num_levels:
        return {}

    # find all rows in the current level with the relevant parent id
    matches = df.loc[df[f"parent_{identifier}"] == parent][identifier].to_list()

    # initialize the dictionary for the current parent
    dict_ = {}

    # for each match, recursively build the dictionary for the next level
    for item in matches:
        dict_[item] = build_chain_dict(df, identifier, item, level + 1, num_levels)

    return dict_


def flatten_dict(d: Dict, parent_keys: list = None, sep: str = "_"):
    """
    Recursively flattens a nested dictionary into a list of dictionaries,
    where each dictionary represents a flattened row of the original dictionary.

    Args:
        d (dict): The input dictionary to be flattened.
        parent_keys (list, optional): The list of parent keys for the current level
            of recursion. Defaults to an empty list.
        sep (str, optional): The separator to be used between parent keys and child
            keys. Defaults to "_".

    Returns:
        list: A list of dictionaries, where each dictionary represents a flattened row
            of the original dictionary.
    """
    if parent_keys is None:
        parent_keys = []
    rows = []
    for k, v in d.items():
        new_keys = parent_keys + [k]
        if isinstance(v, dict) and v:
            rows.extend(flatten_dict(v, new_keys, sep=sep))
        else:
            row = {f"level_{i}": key for i, key in enumerate(new_keys)}
            rows.append(row)
    return rows


def breaks_chain(row, num_levels: int):
    """
    Determines if there is a break in the chain of intents for a given row.
    Args:
      row (pd.Series): A row from a DataFrame containing intent columns.
      num_levels (int): The number of intent levels to consider. If greater than 2,
        an additional intent level is added.
    Returns:
      bool: True if there is a break in the chain of intents, False otherwise.
    """

    intents = ["intent_0", "intent_1"]
    if num_levels > 2:
        intents.append("intent_2")
    for i in range(len(intents) - 1):
        if (
            (pd.isna(row[intents[i]]) or row[intents[i]] == "N/A")
            and not pd.isna(row[intents[i + 1]])
            and row[intents[i + 1]] != "N/A"
        ):
            return True
    return False


def ensure_levels(df: pd.DataFrame, num_levels: int) -> pd.DataFrame:
    """Ensure that the DataFrame has the specified number of levels."""
    for i in range(num_levels):
        if f"level_{i}" not in df.columns:
            df[f"level_{i}"] = np.nan
    return df


def transform_long(
    data: pd.DataFrame,
    intents: list,
    num_levels: int,
) -> pd.DataFrame:
    """
    Transforms the given data by selecting specific levels and intents.

    Args:
        data (pandas.DataFrame): The input data.
        intents (list): The list of intents to filter the data.

    Returns:
        pandas.DataFrame: The transformed data containing selected levels and intents.
    """

    # select all level_0 with an intent_0 that is either result or methodology
    level_0 = data.loc[data["intent_0"].isin(intents)]

    level_0_ids = list(set(level_0["level_0"].to_list()))

    # do the same for level_1 if level_0 is result and methodology
    level_1 = level_0.loc[level_0["intent_1"].isin(intents)]

    level_1_ids = list(set(level_1["level_1"].to_list()))

    if num_levels > 2:
        # do the same for level_2 if level_1 is result and methodology
        level_2 = level_1.loc[level_1["intent_2"].isin(intents)]

        level_2_ids = list(set(level_2["level_2"].to_list()))

    # concatenate the lists of ids
    ids = level_0_ids + level_1_ids
    if num_levels > 2:
        ids += level_2_ids

    # create a list of levels
    levels = [0] * len(level_0_ids) + [1] * len(level_1_ids)

    if num_levels > 2:
        levels += [2] * len(level_2_ids)

    # create a dataframe from the lists of ids and levels
    df_ids = pd.DataFrame({"paper_id": ids, "level": levels})

    df_ids["paper_id"] = df_ids["paper_id"].astype(str)

    df_ids.drop_duplicates(subset=["paper_id"], inplace=True)

    return df_ids


def transform_long_pairs(data, intents, num_levels: int):
    """
    Transforms hierarchical data into a long format DataFrame.

    Args:
      data (pd.DataFrame): The input DataFrame containing hierarchical data.
      intents (list): A list of intents to filter the data at each level.
      num_levels (int): The number of hierarchical levels to process.
    Returns:
      pd.DataFrame: A DataFrame with columns ["paper_id", "parent_paper_id", "intent", "level"].
    """

    df_ids = pd.DataFrame(columns=["paper_id", "parent_paper_id", "intent", "level"])
    data_cp = data.copy()

    for i in range(num_levels):
        level = data_cp.loc[data_cp[f"intent_{i}"].isin(intents)]
        level_ids = level[f"level_{i}"].tolist()
        parent_ids = (
            level[f"level_{str(i-1)}"].tolist() if i > 0 else level["seed_id"].tolist()
        )
        intent = level[f"intent_{i}"].tolist()
        data_cp = level

        df_level = pd.DataFrame(
            {
                "paper_id": level_ids,
                "parent_paper_id": parent_ids,
                "intent": intent,
                "level": [i] * len(level_ids),
            }
        )

        df_ids = pd.concat([df_ids, df_level])

    df_ids["paper_id"] = df_ids["paper_id"].astype(str)
    df_ids.drop_duplicates(
        subset=["parent_paper_id", "paper_id", "level"], inplace=True
    )

    return df_ids
