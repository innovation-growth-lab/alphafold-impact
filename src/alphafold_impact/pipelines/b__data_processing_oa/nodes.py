"""This module contains functions to process data for the OpenAlex pipeline.

The module contains functions to process data by level, fetch additional mesh terms,
and combine data from multiple levels.

Functions:
    - process_data_by_level: Process data by level and return a DataFrame with
        selected columns.
    - fetch_additional_mesh: Fetch additional mesh terms for DOIs in a DataFrame.
    - combine_levels_data: Combine multiple dataframes into a single dataframe and
        remove duplicate rows.
    - process_subfield_data: Process data by level without enriching the data with
        additional mesh terms.
"""

import logging
from typing import Dict, Generator, Tuple
import pandas as pd
import numpy as np
from tqdm import tqdm
from kedro.io import AbstractDataset

logger = logging.getLogger(__name__)


# enriching the data with additional mesh tags
def _data_generator(
    data: Dict[str, AbstractDataset],
    level: int,
) -> Generator[pd.DataFrame, None, None]:
    """
    Load a batch of data from the input dataset.
    """
    # select keys that start by the level
    keys = [key for key in data.keys() if key.startswith(str(level))]

    for key in keys:
        yield _json_loader(data, key, level)


def _json_loader(
    data: Dict[str, AbstractDataset], key: str, level: int
) -> pd.DataFrame:
    """
    Load JSON data, transform it into a DataFrame, and perform data transformations.

    Args:
        data (Dict[str, AbstractDataset]): The input dataset.
        key (str): The key to load from the input dataset.
        level (int): The level of the data.

    Returns:
        pandas.DataFrame: The transformed DataFrame.

    """

    # load json and keep id, doi, publication_date, mesh_terms if exist
    raw_json_data = data[key]()

    output = []

    for parent_id, children_list in raw_json_data.items():

        json_data = [
            {
                k: v
                for k, v in item.items()
                if k
                in [
                    "id",
                    "ids",
                    "doi",
                    "publication_date",
                    "mesh_terms",
                    "cited_by_count",
                    "counts_by_year",
                    "authorships",
                    "topics",
                    "concepts",
                    "fwci",
                    "citation_normalized_percentile",
                    "cited_by_percentile_year",
                ]
            }
            for item in children_list
        ]

        # transform to dataframe and add parent_id
        df = pd.DataFrame(json_data)
        df["parent_id"] = parent_id
        df["level"] = level

        # extract pmid from ids
        df["pmid"] = df["ids"].apply(
            lambda x: (
                x.get("pmid").replace("https://pubmed.ncbi.nlm.nih.gov/", "")
                if x and x.get("pmid")
                else None
            )
        )

        # keep only a list of tuples with "descriptor_ui" and "descriptor_name" for each
        df["mesh_terms"] = df["mesh_terms"].apply(
            lambda x: (
                [(c["descriptor_ui"], c["descriptor_name"]) for c in x] if x else None
            )
        )

        # break atuhorship nested dictionary jsons, create triplets of authorship
        df["authorships"] = df["authorships"].apply(
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

        # change doi to remove the url
        df["doi"] = df["doi"].str.replace("https://doi.org/", "")

        # create a list of topics
        df["topics"] = df["topics"].apply(
            lambda x: (
                [
                    (
                        topic["id"].replace("https://openalex.org/", ""),
                        topic["display_name"],
                        topic["subfield"]["id"].replace("https://openalex.org/", ""),
                        topic["subfield"]["display_name"],
                        topic["field"]["id"].replace("https://openalex.org/", ""),
                        topic["field"]["display_name"],
                        topic["domain"]["id"].replace("https://openalex.org/", ""),
                        topic["domain"]["display_name"],
                    )
                    for topic in x
                ]
                if x
                else None
            )
        )

        # extract concepts
        df["concepts"] = df["concepts"].apply(
            lambda x: (
                [
                    (
                        concept["id"].replace("https://openalex.org/", ""),
                        concept["display_name"],
                    )
                    for concept in x
                ]
                if x
                else None
            )
        )

        # Extract the content of citation_normalized_percentile
        try:
            df[
                [
                    "citation_normalized_percentile_value",
                    "citation_normalized_percentile_is_in_top_1_percent",
                    "citation_normalized_percentile_is_in_top_10_percent",
                ]
            ] = df.apply(
                lambda x: (pd.Series(x["citation_normalized_percentile"])),
                axis=1,
                result_type="expand",
            )
        except ValueError:
            logger.warning(
                "citation_normalized_percentile not found in %s", df["id"].values[0]
            )

        # Extract the content of cited_by_percentile_year
        df[
            [
                "cited_by_percentile_year_min",
                "cited_by_percentile_year_max",
            ]
        ] = df.apply(
            lambda x: pd.Series(x["cited_by_percentile_year"]),
            axis=1,
            result_type="expand",
        )

        # remove any column that is all NAN
        df.dropna(axis=1, how="all", inplace=True)

        # append to output
        output.append(df)

    df = pd.concat(output, ignore_index=True)

    return df


def process_data_by_level(data: Dict[str, AbstractDataset], level: int) -> pd.DataFrame:
    """
    Process data by level and return a DataFrame with selected columns.

    Args:
        data (Dict[str, AbstractDataset]): A dictionary containing the input data.
        level (int): The level to process the data for.
        extra_mesh (bool): Whether to enrich the data with additional mesh terms.

    Returns:
        pd.DataFrame: A DataFrame containing the processed data with selected columns.
    """

    logger.info("Processing data for level %s", level)

    # Generate data for current level
    data_gen = _data_generator(data, level)

    # Initialize an empty DataFrame for current level
    df_level = pd.DataFrame()

    # Iterate over data batches in current level
    for df_batch in data_gen:

        # Append current batch to level DataFrame
        df_level = pd.concat([df_level, df_batch])

    return df_level[
        [
            "parent_id",
            "id",
            "pmid",
            "level",
            "doi",
            "publication_date",
            "mesh_terms",
            "cited_by_count",
            "counts_by_year",
            "authorships",
            "topics",
            "concepts",
            "fwci",
            "citation_normalized_percentile_value",
        ]
    ]


def combine_levels_data(unique: str = "all", **kwargs) -> pd.DataFrame:
    """
    Combines multiple dataframes into a single dataframe and removes duplicate rows
        based on the 'id' column.

    Parameters:
    unique (str): The method to remove duplicate rows. Options are 'id', 'level', or 'all'.
    **kwargs: keyword arguments representing the dataframes to be combined.

    Returns:
    pd.DataFrame: A combined dataframe with duplicate rows removed based
        on the 'id' column.
    """
    assert unique in [
        "id",
        "level",
        "all",
    ], f"unique must be either 'id' or 'level, or 'all'. Instead, got {unique}."
    df = pd.concat([level for level in kwargs.values()])
    if unique == "id":
        logger.info("Removing duplicate rows based on 'id' column.")
        df = df.sort_values(["level", "id"], ascending=True).drop_duplicates(
            ["id"], keep="first"
        )
    elif unique == "level":
        logger.info("Keeping duplicate rows based on 'id' column.")
        # keep only the duplicated ids that are the lowest level, but keep all of those
        df = df.sort_values(["level", "id"], ascending=True)
        lowest_levels = df.groupby("id")["level"].min()
        df = df.merge(lowest_levels, on="id", suffixes=["", "_min"])
        df = df[df["level"] == df["level_min"]]
        df = df.drop(columns="level_min")
    elif unique == "all":
        logger.info("Removing only full duplicate rows.")
        df = df.drop_duplicates(
            subset=[
                col
                for col in df.columns
                if col
                not in [
                    "strength",
                    "mesh_terms",
                    "authorships",
                    "topics",
                    "concepts",
                    "counts_by_year",
                ]
            ]
        )

        df.reset_index(drop=True, inplace=True)
        df["level"] = df["level"].astype(str)

    return df


def combine_levels_data_counterfactuals(unique: str = "all", **kwargs) -> pd.DataFrame:
    """
    Combines multiple dataframes into a single dataframe and removes fully duplicate rows.

    Parameters:
    unique (str): The method to remove duplicate rows. Options are 'id', 'level', or 'all'.
    **kwargs: keyword arguments representing the dataframes to be combined.

    Returns:
    pd.DataFrame: A combined dataframe with duplicate rows removed based
        on the 'id' column.
    """
    assert unique in [
        "id",
        "level",
        "all",
    ], f"unique must be either 'id' or 'level, or 'all'. Instead, got {unique}."

    # filter out rows in level t+1 if their parent_id is not in level t
    levels = [level for level in kwargs.values()]
    for i in range(len(levels) - 1):
        levels[i + 1] = levels[i + 1][levels[i + 1]["parent_id"].isin(levels[i]["id"])]

    df = pd.concat([level for level in levels])

    logger.info("Removing only full duplicate rows.")
    df = df.drop_duplicates(
        subset=[
            col
            for col in df.columns
            if col
            not in [
                "strength",
                "mesh_terms",
                "authorships",
                "topics",
                "concepts",
                "counts_by_year",
            ]
        ]
    )

    df.reset_index(drop=True, inplace=True)
    df["level"] = df["level"].astype(str)

    return df


def process_data_by_level_ptd(
    data: Dict[str, AbstractDataset], level: int
) -> Generator[Dict[str, pd.DataFrame], None, None]:
    """
    Process data by level and return a DataFrame with selected columns.

    Args:
        data (Dict[str, AbstractDataset]): A dictionary containing the input data.
        level (int): The level to process the data for.
        extra_mesh (bool): Whether to enrich the data with additional mesh terms.

    Returns:
        pd.DataFrame: A DataFrame containing the processed data with selected columns.
    """

    logger.info("Processing data for level %s", level)

    # Generate data for current level
    data_gen = _data_generator(data, level)

    for i, df_batch in enumerate(
        tqdm(data_gen, total=len(data), desc="Processing data partitions")
    ):
        logger.info("Processing data partition: %s / %s", i + 1, len(data))
        yield {
            f"s{i}": df_batch[
                [
                    "parent_id",
                    "id",
                    "pmid",
                    "level",
                    "doi",
                    "publication_date",
                    "mesh_terms",
                    "cited_by_count",
                    "counts_by_year",
                    "authorships",
                    "topics",
                    "concepts",
                    "fwci",
                    "citation_normalized_percentile_value",
                ]
            ]
        }


def concat_pq_ptd(
    data: Dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """
    Concatenate dataframes from multiple levels into a single dataframe.

    Args:
        data (Dict[str, pd.DataFrame]): A dictionary containing the input data.

    Returns:
        pd.DataFrame: A DataFrame containing the concatenated data.
    """
    output = []
    for i, loader in enumerate(data.values()):
        logger.info("Processing data partition: %s / %s", i + 1, len(data))
        data_pt = loader()
        # drop rows with id or parent id W3177828909, W3211795435, W3202105508
        data_pt = data_pt[
            ~(
                data_pt["id"].isin(["W3177828909", "W3211795435", "W3202105508"])
                | data_pt["parent_id"].isin(
                    ["W3177828909", "W3211795435", "W3202105508"]
                )
            )
        ]
        output.append(data_pt)
    return pd.concat(output)


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
                if k
                in [
                    "id",
                    "ids",
                    "doi",
                    "publication_date",
                    "mesh_terms",
                    "cited_by_count",
                    "counts_by_year",
                    "authorships",
                    "topics",
                    "concepts",
                    "fwci",
                    "citation_normalized_percentile",
                    "cited_by_percentile_year",
                ]
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

        # extract pmid from ids
        df["pmid"] = df["ids"].apply(
            lambda x: (
                x.get("pmid").replace("https://pubmed.ncbi.nlm.nih.gov/", "")
                if x and x.get("pmid")
                else None
            )
        )

        # break atuhorship nested dictionary jsons, create triplets of authorship
        df["authorships"] = df["authorships"].apply(
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

        # create a list of topics
        df["topics"] = df["topics"].apply(
            lambda x: (
                [
                    (
                        topic["id"].replace("https://openalex.org/", ""),
                        topic["display_name"],
                        topic["subfield"]["id"].replace("https://openalex.org/", ""),
                        topic["subfield"]["display_name"],
                        topic["field"]["id"].replace("https://openalex.org/", ""),
                        topic["field"]["display_name"],
                        topic["domain"]["id"].replace("https://openalex.org/", ""),
                        topic["domain"]["display_name"],
                    )
                    for topic in x
                ]
                if x
                else None
            )
        )

        # extract concepts
        df["concepts"] = df["concepts"].apply(
            lambda x: (
                [
                    (
                        concept["id"].replace("https://openalex.org/", ""),
                        concept["display_name"],
                    )
                    for concept in x
                ]
                if x
                else None
            )
        )

        # Extract the content of citation_normalized_percentile
        try:
            df[
                [
                    "citation_normalized_percentile_value",
                    "citation_normalized_percentile_is_in_top_1_percent",
                    "citation_normalized_percentile_is_in_top_10_percent",
                ]
            ] = df.apply(
                lambda x: (pd.Series(x["citation_normalized_percentile"])),
                axis=1,
                result_type="expand",
            )
        except ValueError:
            logger.warning(
                "citation_normalized_percentile not found in %s", df["id"].values[0]
            )

        # Extract the content of cited_by_percentile_year
        df[
            [
                "cited_by_percentile_year_min",
                "cited_by_percentile_year_max",
            ]
        ] = df.apply(
            lambda x: pd.Series(x["cited_by_percentile_year"]),
            axis=1,
            result_type="expand",
        )

        # change doi to remove the url
        df["doi"] = df["doi"].str.replace("https://doi.org/", "")

        # change id to remove openalex url
        df["id"] = df["id"].str.replace("https://openalex.org/", "")

        # Append current batch to level DataFrame
        subfield_df = pd.concat([subfield_df, df])

    return subfield_df[
        [
            "id",
            "pmid",
            "doi",
            "publication_date",
            "mesh_terms",
            "cited_by_count",
            "counts_by_year",
            "authorships",
            "concepts",
            "topics",
            "fwci",
            "citation_normalized_percentile_value",
        ]
    ]


def reassign_ct_levels(
    data: pd.DataFrame,
    ct_data: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reassigns the ct levels in the given data based on the provided ct_data.

    Args:
        data (pd.DataFrame): The input data containing the levels to be reassigned.
        ct_data (pd.DataFrame): The ct data used for reassigning the levels.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two dataframes:
            - data_ct: The data with reassigned ct levels.
            - data_other: The remaining data with reassigned levels.

    """
    # drop any rows with id or parent_id W3177828909, W3211795435, W3202105508
    data = data[
        ~(
            data["id"].isin(["W3177828909", "W3211795435", "W3202105508"])
            | data["parent_id"].isin(["W3177828909", "W3211795435", "W3202105508"])
        )
    ]

    data_cp = data.copy()

    logger.info("Seed - Identify the ct_data ids that appear as ids")
    data_cp["ct_seed"] = data_cp.apply(
        lambda row: row["id"] in ct_data["parent_id"].to_list()
        and ((pd.isna(row["level"])) | (row["level"] == "nan")),
        axis=1,
    )

    logger.info("Level 0 - Identify the ct_data ids that appear as parent_id")
    data_cp["ct_l0"] = np.where(
        (data_cp["level"] == "0.0") & (data_cp["parent_id"].isin(ct_data["parent_id"])),
        True,
        False,
    )

    logger.info("Level 1 - Identify the level 1 ids associated to a CT")
    data_cp["ct_l1"] = np.where(
        (data_cp["level"] == "1.0")
        & (data_cp["parent_id"].isin(data_cp[data_cp["ct_l0"]]["id"])),
        True,
        False,
    )

    logger.info("Reassigning the ct levels")
    data_ct = data_cp[data_cp["ct_seed"] | data_cp["ct_l0"] | data_cp["ct_l1"]]
    data_other = data_cp[~(data_cp["ct_seed"] | data_cp["ct_l0"] | data_cp["ct_l1"])]

    logger.info("Reassign levels in CT (seed, 0, 1)")
    level_mapping = {"nan": "seed", "0.0": "0", "1.0": "1"}
    data_ct["level"] = data_ct["level"].map(level_mapping).fillna(data_ct["level"])

    logger.info("Reassign levels in Other (0, 1, 2)")
    level_mapping = {"nan": "0", "0.0": "1", "1.0": "2"}
    data_other["level"] = (
        data_other["level"].map(level_mapping).fillna(data_other["level"])
    )

    # drop columns
    data_ct = data_ct.drop(columns=["ct_seed", "ct_l0", "ct_l1"])
    data_other = data_other.drop(columns=["ct_seed", "ct_l0", "ct_l1"])

    # filter out publications older than 2015 for data_other level 0
    level_0_data = data_other[
        (data_other["level"] == "0") & (data_other["publication_date"] >= "2015-01-01")
    ]

    # HACK: keep 50% of the level_0_data for other SB papers
    level_0_data = level_0_data.sample(frac=0.5, random_state=42)

    # Get level 1 data where parent_id is in level 0 ids
    level_1_data = data_other[
        (data_other["level"] == "1") & data_other["parent_id"].isin(level_0_data["id"])
    ]

    # Get level 2 data
    level_2_data = data_other[
        (data_other["level"] == "2") & data_other["parent_id"].isin(level_1_data["id"])
    ]

    # Concatenate the data
    data_other = pd.concat([level_0_data, level_1_data, level_2_data])

    return data_ct, data_other
