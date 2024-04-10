"""This module contains functions to process iCite data.

Functions:
    - filter_and_combine_icite: Loads iCite data from a partitioned dataset,
        filters out records before a specified year and selects only relevant columns.
"""

from typing import Dict, Callable, Any
import logging
import pandas as pd

logger = logging.getLogger(__name__)


def filter_and_combine_icite(
    icite_partitions: Dict[str, Callable[[], Any]], year=2016
) -> pd.DataFrame:
    """Loads iCite data from a partitioned dataset, filters out
    records before a specified year and selects only relevant columns.

    Args:
        icite_partitions (dict): iCite partitioned dataset.
        year (int): Year to filter out records before.

    Returns:
        DataFrame of combined and filtered iCite data.
    """
    filtered_icite = []

    for name, icite_partition in icite_partitions.items():
        logger.info("Filtering iCite partition: %s", name)
        icite_partition = icite_partition()
        current_matches = icite_partition.query(f"year >= {year}")
        logger.info("Adding %s matches", len(current_matches))
        filtered_icite.append(current_matches)

    return pd.concat(filtered_icite)[
        [
            "pmid",
            "doi",
            "title",
            "authors",
            "year",
            "journal",
            "is_research_article",
            "citation_count",
            "apt",
            "is_clinical",
            "cited_by_clin",
            "cited_by",
            "references",
        ]
    ].reset_index(drop=True)
