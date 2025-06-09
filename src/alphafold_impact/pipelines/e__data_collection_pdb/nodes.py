"""
This is a boilerplate pipeline 'data_collection_pdb'
generated using Kedro 0.19.1
"""

import logging
import asyncio
from typing import Dict
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from .utils import (
    fetch_pdb_ids,
    create_pdb_query,
    collect_structure_matches_async,
)

logger = logging.getLogger(__name__)


def fetch_pbd_details(config: Dict[str, str]) -> pd.DataFrame:
    """Fetch the details of a PDB id"""
    pbd_ids = fetch_pdb_ids()
    session = requests.Session()
    retries = Retry(
        total=config["max_retries"], backoff_factor=config["backoff_factor"]
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))
    url = "https://data.rcsb.org/graphql"

    # create sublists of 50 items
    pbd_ids = [pbd_ids[i : i + 50] for i in range(0, len(pbd_ids), 50)]

    outputs = []
    for ids in pbd_ids:
        logger.info("Fetching details for %s", ids)
        query = create_pdb_query(ids)
        response = requests.get(url, params={"query": query}, timeout=10)
        data = response.json()
        data = pd.DataFrame(data["data"]["entries"])

        # get date
        data["rcsb_accession_info"] = data["rcsb_accession_info"].apply(
            lambda x: x["initial_release_date"]
        )

        # get an author list from the list of dictionaries with key "name"
        data["audit_authors"] = data["audit_author"].apply(
            lambda x: ", ".join([d.get("name", "") for d in x])
        )

        # get pmid
        data["pmid"] = data["rcsb_primary_citation"].apply(
            lambda x: x.get("pdbx_database_id_PubMed", "") if x is not None else ""
        )

        # get doi
        data["doi"] = data["rcsb_primary_citation"].apply(
            lambda x: x.get("pdbx_database_id_DOI", "") if x is not None else ""
        )

        # extract resolution
        data["resolution"] = data["rcsb_entry_info"].apply(
            lambda x: (
                x.get("resolution_combined", [""])[0]
                if x.get("resolution_combined")
                else ""
            )
        )
        # extract R-free factor
        data["R_free"] = data["refine"].apply(
            lambda x: (
                x[0].get("ls_R_factor_R_free", "") if isinstance(x, list) else None
            )
        )

        # append
        outputs.append(data)

    # concatenate
    outputs = pd.concat(outputs)

    # make sure ids are strings
    outputs["rcsb_id"] = outputs["rcsb_id"].astype(str)
    outputs["pmid"] = outputs["pmid"].astype(str)
    outputs["doi"] = outputs["doi"].astype(str)
    outputs["resolution"] = outputs["resolution"].astype(str)
    outputs["R_free"] = outputs["R_free"].astype(str)

    return outputs


def collect_structure_match_scores(
    pdb_entries: pd.DataFrame, config: Dict[str, str]
) -> pd.DataFrame:
    """
    Collect structure match scores for PDB entries using the RCSB PDB search API.

    This function queries the PDB structure search API asynchronously for each entry
    to find similar structures and their match scores. It returns the top 10 matches
    for each query structure with proper rate limiting and retry logic.

    Args:
        pdb_entries (pd.DataFrame): DataFrame containing PDB entries with 'rcsb_id' column.
        config (Dict[str, str]): Configuration parameters including:
            - calls_per_second: Rate limit for API calls (default: 2.0)
            - max_concurrent: Maximum concurrent requests (default: 10)
            - max_retries: Maximum retries per request (default: 5)

    Returns:
        pd.DataFrame: DataFrame with columns 'query', 'target', 'score' containing
                     structure match information.
    """
    logger.info(
        "Starting async structure match score collection for %d PDB entries",
        len(pdb_entries),
    )

    # Get list of unique PDB IDs to query
    pdb_ids = pdb_entries["rcsb_id"].unique().tolist()

    # Configuration parameters
    calls_per_second = float(config.get("calls_per_second", 2))
    max_concurrent = int(config.get("max_concurrent", 10))
    max_retries = int(config.get("max_retries", 5))

    logger.info(
        "Configuration: %.1f calls/sec, %d max concurrent, %d max retries",
        calls_per_second,
        max_concurrent,
        max_retries,
    )

    # Run the async collection
    result_df = asyncio.run(
        collect_structure_matches_async(
            pdb_ids=pdb_ids,
            calls_per_second=calls_per_second,
            max_concurrent=max_concurrent,
            max_retries=max_retries,
        )
    )

    return result_df
