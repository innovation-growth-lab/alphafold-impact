"""
This module contains functions for processing funding opportunities and
related MeSH data.
"""

import pandas as pd


def _find_lower_lvl_mesh_terms(
    full_mesh_tree: pd.DataFrame, tree_numbers: list
) -> pd.DataFrame:
    """
    Filter the DataFrame by rows where the 'tree_number' column starts with any of
    the strings given in tree_numbers.

    Args:
        full_mesh_tree: A pandas DataFrame with a 'tree_number' column.
        tree_number: A string to match the start of the 'tree_number' entries.

    Returns:
        A DataFrame containing the rows where 'tree_number' starts with the given string.
    """
    return full_mesh_tree[
        full_mesh_tree["tree_number"].apply(
            lambda x: any(
                x.startswith(tree_number) for tree_number in tree_numbers if pd.notna(x)
            )
        )
    ]


def _list_mesh_duis(subfield_mesh_dict: dict, full_mesh_tree: pd.DataFrame) -> list:
    """
    Return list of MeSH DUIs

    Args:
        subfield_mesh_dict (dict): Dictionary with keys of MeSH names and
            list of MeSH tree numbers for the values
        full_mesh_tree (dict): Full collection of MeSH names,
            DUIs and tree numbers
    """
    return [
        DUI
        for tree_numbers in subfield_mesh_dict.values()
        for DUI in _find_lower_lvl_mesh_terms(full_mesh_tree, tree_numbers).DUI.values
    ]


def _search_text_in_multiple_cols(
    df: pd.DataFrame, columns: list[str], search_terms: list[str]
) -> pd.DataFrame:
    """
    Search for any of several texts in any of multiple specified columns of
    a pandas DataFrame and return rows containing those texts, ignoring case
    differences.

    Args:
        df (pd.DataFrame): The DataFrame to search in.
        columns (list[str]): The names of the columns to perform the search.
        texts (list[str]): The texts to search for.

    Returns:
        pd.DataFrame: A DataFrame containing the rows where any text is found
            in any column, case-insensitively.
    """
    matching_rows = pd.DataFrame()
    for column in columns:
        if column in df.columns:
            for search_term in search_terms:
                current_matches = df[
                    df[column].str.contains(search_term, case=False, na=False)
                ]
                matching_rows = pd.concat(
                    [matching_rows, current_matches]
                ).drop_duplicates()
    return matching_rows


def find_relevant_funding_opps(
    nih_funding_opps: pd.DataFrame,
    nih_funding_opps_mesh_tags: pd.DataFrame,
    full_mesh_tree: pd.DataFrame,
    subfield_label: str,
    subfield_mesh_dict: dict,
    search_terms: list,
    search_cols: list,
) -> pd.DataFrame:
    # Find relevant funding opps based on MeSH term
    mesh_duis = _list_mesh_duis(subfield_mesh_dict, full_mesh_tree)
    nih_funding_opps_mesh_tags["mesh_term"] = nih_funding_opps_mesh_tags[
        "mesh_term"
    ].str.replace("*", "", regex=False)
    relevant_nih_funding_opps_mesh_tags = nih_funding_opps_mesh_tags.query(
        f"mesh_term_id in {mesh_duis}"
    )
    nih_funding_opp_ids_from_mesh = [
        int(text_id)
        for text_id in relevant_nih_funding_opps_mesh_tags.text_id.to_list()
    ]
    nih_funding_opps_from_mesh = nih_funding_opps.query(
        f"id in {nih_funding_opp_ids_from_mesh}"
    ).assign(source="nih mesh", subfield=subfield_label)
    # Find relevant funding opps based on keyword search
    nih_funding_opps_from_search = _search_text_in_multiple_cols(
        nih_funding_opps, search_cols, search_terms
    ).assign(source="nih search", subfield=subfield_label)
    # Combine funding opps found from MeSH terms and keyword search
    return (
        pd.concat([nih_funding_opps_from_mesh, nih_funding_opps_from_search])
        .drop_duplicates(subset="id")
        .reset_index(drop=True)
    )
