"""This module contains functions to process iCite data.

Functions:
    - filter_and_combine_icite: Loads iCite data from a partitioned JSON dataset,
        filters out records before a specified year and selects only relevant columns.
"""

from typing import Dict, Callable, Any
import logging
import pandas as pd

logger = logging.getLogger(__name__)


def filter_and_combine_icite(
    icite_partitions: Dict[str, Callable[[], Any]], year=2015
) -> pd.DataFrame:
    """Loads iCite data from a partitioned JSON dataset, filters out
    records before a specified year and selects only relevant columns.

    Args:
        icite_partitions (dict): iCite partitioned JSON dataset.
        year (int): Year to filter out records before.

    Returns:
        DataFrame of combined and filtered iCite data.
    """
    filtered_icite = []
    desired_columns = [
        "pmid",
        "doi",
        "cited_by_clin",
    ]

    for name, icite_partition in icite_partitions.items():
        logger.info("Processing iCite partition: %s", name)

        # Load the JSON data from the partition
        json_data = icite_partition()

        # Convert JSON data to DataFrame
        if isinstance(json_data, list):
            # If it's a list of records
            icite_df = pd.DataFrame(json_data)
        elif isinstance(json_data, dict):
            # If it's a single record or nested structure
            icite_df = pd.json_normalize(json_data)
        else:
            logger.warning("Unexpected data format in partition %s", name)
            continue

        if icite_df.empty:
            logger.warning("No data found in partition %s", name)
            continue

        # Filter by year if the year column exists
        if "year" in icite_df.columns:
            current_matches = icite_df.query(f"year >= {year}")
            logger.info(
                "Filtered %d records from partition %s (year >= %d)",
                len(current_matches),
                name,
                year,
            )
        else:
            logger.warning(
                "Year column not found in partition %s, including all records", name
            )
            current_matches = icite_df

        if not current_matches.empty:
            available_columns = [
                col for col in desired_columns if col in current_matches.columns
            ]
            filtered_icite.append(current_matches[available_columns])

    if not filtered_icite:
        logger.warning("No data found in any partitions")
        return pd.DataFrame()

    # Combine all partitions
    result_df = pd.concat(filtered_icite, ignore_index=True)

    return result_df
