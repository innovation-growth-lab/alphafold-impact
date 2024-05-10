"""
This is a boilerplate pipeline 'data_processing_chains'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict
import pandas as pd
import numpy as np
from Bio import Entrez

pd.options.mode.copy_on_write = True

Entrez.email = "david.ampudia@nesta.org.uk"
logger = logging.getLogger(__name__)


def _sort_and_drop(
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

    # data = data[~data["intent"].isin(["", None])]

    # keep the first row of every parent_{identifier}, identifier
    data = data.groupby(["parent_id", "id"]).first().reset_index()

    if unique:
        logger.info("Dropping duplicate citation links")
        data.sort_values("sort_order").drop_duplicates(subset="id", inplace=True)

    # drop the sort_order column, reindex
    data.drop(columns="sort_order", inplace=True)
    data.reset_index(drop=True, inplace=True)

    return data


def _build_chain_dict(df, identifier, parent, level):
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
    if level > 2:
        return {}

    # find all rows in the current level with the relevant parent id
    matches = df.loc[df[f"parent_{identifier}"] == parent][identifier].to_list()

    # initialize the dictionary for the current parent
    dict_ = {}

    # for each match, recursively build the dictionary for the next level
    for item in matches:
        dict_[item] = _build_chain_dict(df, identifier, item, level + 1)

    return dict_


def _flatten_dict(d: Dict, parent_keys: list = None, sep: str = "_"):
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
            rows.extend(_flatten_dict(v, new_keys, sep=sep))
        else:
            row = {f"level_{i}": key for i, key in enumerate(new_keys)}
            rows.append(row)
    return rows


def _breaks_chain(row):
    intents = ["intent_0", "intent_1", "intent_2"]
    for i in range(len(intents) - 1):
        if (
            (pd.isna(row[intents[i]]) or row[intents[i]] == "N/A")
            and not pd.isna(row[intents[i + 1]])
            and row[intents[i + 1]] != "N/A"
        ):
            return True
    return False


def filter_relevant_citation_links(
    alphafold_data: pd.DataFrame,
    identifier: str,
    **kwargs,
) -> pd.DataFrame:
    """
    Filters and processes citation links based on relevance.

    Args:
        alphafold_data (pd.DataFrame): The input DataFrame containing citation
             link data.
        identifier (str): The identifier used to filter the citation links.
        **kwargs: Additional keyword arguments for custom sorting.

    Returns:
        pd.DataFrame: The processed DataFrame containing relevant citation chains.
    """

    logger.info("Exploding citation links")
    alphafold_data = alphafold_data.explode("strength")

    # logger.info("Remove citation links with no strength")
    # alphafold_data.dropna(subset=["strength"], inplace=True)

    for col in ["intent", "context"]:
        logger.info("Extracting %s from citation links", col)
        alphafold_data[col] = alphafold_data["strength"].apply(
            lambda x: x[col] if x else None  # pylint: disable=cell-var-from-loop
        )

    # create strong variable
    alphafold_data["intent"] = alphafold_data["intent"].apply(
        lambda x: (
            "strong"
            if x in ["result", "methodology"]
            else (
                "weak" if x in ["background"] else "unknown" if x in [""] else "no_data"
            )
        )
    )

    # define & implement the custom sorting order
    alphafold_data = _sort_and_drop(alphafold_data, **kwargs)

    if identifier == "pmid":
        # assign 99999999 PMID to Multimer paper, with OA id W3202105508
        alphafold_data.loc[
            alphafold_data["parent_id"] == "W3202105508", "parent_pmid"
        ] = "99999999"

    logger.info("Filter nulls for relevant identifier %s", identifier)
    alphafold_data = alphafold_data[alphafold_data[identifier].notnull()]

    logger.info("Drop duplicates of parent_id, identifier, level, and intent")
    alphafold_data.drop_duplicates(
        subset=["parent_id", identifier, "intent"], inplace=True
    )

    output = []
    for seed_id in list(
        alphafold_data[alphafold_data["level"] == 0][f"parent_{identifier}"].unique()
    ):
        logger.info("Building citation chain for seed %s", seed_id)
        data_dict = _build_chain_dict(alphafold_data, identifier, seed_id, 0)
        flat_rows = _flatten_dict(data_dict)
        seed_df = pd.DataFrame(flat_rows)
        seed_df["seed_id"] = seed_id
        seed_df = seed_df[
            [
                "seed_id",
                "level_0",
                "level_1",
                "level_2",
            ]
        ]
        output.append(seed_df)

    chains_df = pd.concat(output)

    # for each level, merge chains with paper intent data
    for i in range(4):
        if i == 0:
            chains_df = chains_df.merge(
                alphafold_data[[f"parent_{identifier}", identifier, "intent"]],
                left_on=["seed_id", "level_0"],
                right_on=[f"parent_{identifier}", identifier],
                how="left",
            )
        else:
            chains_df = chains_df.merge(
                alphafold_data[[f"parent_{identifier}", identifier, "intent"]],
                left_on=[f"level_{i-1}", f"level_{i}"],
                right_on=[f"parent_{identifier}", identifier],
                how="left",
            )
        # rename the intent column
        chains_df.rename(columns={"intent": f"intent_{i}"}, inplace=True)

        # drop duplicate columns
        chains_df.drop(columns=[f"parent_{identifier}", identifier], inplace=True)

    # substitute empty string in intent for None
    for i in range(4):
        chains_df[f"intent_{i}"] = chains_df[f"intent_{i}"].apply(
            lambda x: x if x != "" else None
        )

    # in the intent columns, replace NaN with 'N/A'
    for i in range(4):
        chains_df[f"intent_{i}"] = chains_df[f"intent_{i}"].apply(
            lambda x: "N/A" if x is np.nan else x
        )
    complete_chains_df = chains_df[~chains_df.apply(_breaks_chain, axis=1)]

    # drop rows that have all four intent as nan or N/A
    complete_chains_df = complete_chains_df[
        ~complete_chains_df[["intent_0", "intent_1", "intent_2"]]
        .applymap(lambda x: x == "N/A" or pd.isna(x))
        .all(axis=1)
    ]

    return complete_chains_df


def _transform_long(data: pd.DataFrame, intents: list) -> pd.DataFrame:
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

    # do the same for level_2 if level_1 is result and methodology
    level_2 = level_1.loc[level_1["intent_2"].isin(intents)]

    level_2_ids = list(set(level_2["level_2"].to_list()))

    # concatenate the lists of ids
    ids = level_0_ids + level_1_ids + level_2_ids

    # create a list of levels
    levels = [0] * len(level_0_ids) + [1] * len(level_1_ids) + [2] * len(level_2_ids)

    # create a dataframe from the lists of ids and levels
    df_ids = pd.DataFrame({"paper_id": ids, "level": levels})

    df_ids["paper_id"] = df_ids["paper_id"].astype(str)

    df_ids.drop_duplicates(subset=["paper_id"], inplace=True)

    return df_ids


def get_entrez_ptype_pmid(pmid: str) -> str:
    """
    Retrieves the publication type for a given PubMed ID (PMID).

    Args:
        pmid (str): The PubMed ID (PMID) of the publication.

    Returns:
        str: The publication type of the specified publication.
    """
    stream = Entrez.efetch(db="pubmed", id=pmid, retmax="1")
    record = Entrez.read(stream)
    return str(
        record["PubmedArticle"][0]["MedlineCitation"]["Article"].get(
            "PublicationTypeList"
        )[0]
    )


def get_papers_with_clinical_article_citations(
    chains: pd.DataFrame,
    icite_data: pd.DataFrame,
    identifier: str,
) -> pd.DataFrame:
    """
    Retrieves papers with clinical article citations based on the given identifier.

    Args:
        chains (pd.DataFrame): The chains dataframe.
        icite_data (pd.DataFrame): The iCite data dataframe.
        identifier (str): The identifier to use for matching (e.g., "pmid").

    Returns:
        pd.DataFrame: The dataframe containing papers with clinical article citations.
    """
    if identifier == "pmid":
        icite_data["pmid"] = icite_data["pmid"].astype(str)

    logger.info("Transforming long chains for identifier %s", identifier)
    data_ids = _transform_long(chains, ["strong"])

    logger.info("Merging chains with iCite data")
    data_ids["cited_by_clin"] = data_ids["paper_id"].apply(
        lambda x: (
            icite_data.loc[icite_data[identifier] == x]["cited_by_clin"].values[0]
            if x in icite_data[identifier].to_list()
            else None
        )
    )

    logger.info("Exploding cited_by_clin")
    data_ids = data_ids[data_ids["cited_by_clin"].astype(str) != "nan"]

    # drop if cited_by_clin is None
    data_ids = data_ids[data_ids["cited_by_clin"].notnull()]

    data_ids["cited_by_clin"] = data_ids["cited_by_clin"].apply(lambda x: x.split(" "))
    data_ids = data_ids.explode("cited_by_clin")

    logger.info("Creating clinical article links")
    data_ids["ca_link"] = data_ids["cited_by_clin"].apply(
        lambda x: f"https://pubmed.ncbi.nlm.nih.gov/{x}"
    )

    # reindex
    data_ids.reset_index(inplace=True, drop=True)

    # drop dup
    data_ids.drop_duplicates(subset=["paper_id", "cited_by_clin"], inplace=True)

    if identifier == "pmid":
        data_ids["publication_type"] = data_ids["cited_by_clin"].apply(
            get_entrez_ptype_pmid
        )

    return data_ids


def _transform_long_pairs(data, intents):
    df_ids = pd.DataFrame(columns=["paper_id", "parent_paper_id", "intent", "level"])
    data_cp = data.copy()

    for i in range(4):
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


def get_papers_with_full_chain(
    chains: pd.DataFrame,
    alphafold_data: pd.DataFrame,
    identifier: str,
) -> pd.DataFrame:
    """
    Retrieves papers with the full chain of data for a given identifier.

    Args:
        chains (pd.DataFrame): DataFrame containing the chains data.
        alphafold_data (pd.DataFrame): DataFrame containing the AlphaFold data.
        identifier (str): Identifier used to match the chains and AlphaFold data.

    Returns:
        pd.DataFrame: DataFrame containing the merged data of chains and AlphaFold data.
    """
    logger.info("Transforming long chains for identifier %s", identifier)
    data_ids = _transform_long(chains, ["strong"])

    logger.info("Merging chains with AlphaFold data")
    alphafold_data = alphafold_data.merge(
        data_ids, how="inner", left_on=identifier, right_on="paper_id"
    )
    return alphafold_data


def get_chain_label_papers(
    chains: pd.DataFrame,
    alphafold_data: pd.DataFrame,
    identifier: str,
) -> pd.DataFrame:
    """
    Identifies and labels chains in the alphafold_data DataFrame based on the information
    in the chains DataFrame.

    Args:
        chains (pd.DataFrame): DataFrame containing information about chains.
        alphafold_data (pd.DataFrame): DataFrame containing alphafold data.
        identifier (str): Identifier column name in the alphafold_data DataFrame.

    Returns:
        pd.DataFrame: DataFrame with the labeled chains.

    Raises:
        None
    """
    logger.info("Identifying strong-only chains")
    data_ids_strong = _transform_long_pairs(chains, ["strong"])

    logger.info("Identifying weak-only chains")
    data_ids_weak = _transform_long_pairs(chains, ["weak"])

    logger.info("Identifying unknown only chains")
    data_ids_unknown = _transform_long_pairs(chains, ["unknown"])

    logger.info("Identifying no_data chains only")
    data_ids_no_data = _transform_long_pairs(chains, ["no_data"])

    logger.info("Identifying partial only-strong chains")
    data_ids_partial_strong = _transform_long_pairs(
        chains, ["strong", "no_data", "unknown"]
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
    data_ids_partial_weak = _transform_long_pairs(
        chains, ["weak", "no_data", "unknown"]
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
    data_ids_mixed = _transform_long_pairs(chains, ["strong", "weak"])

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

    # map chain_label to alphafold_data identifier, parent_{identifier}
    alphafold_data["chain_label"] = alphafold_data.apply(
        lambda x: chain_labels_map.get(
            (x[identifier], x[f"parent_{identifier}"]), None
        ),
        axis=1,
    )

    # fillna chain_label with "other"
    alphafold_data["chain_label"].fillna("other", inplace=True)

    return alphafold_data
