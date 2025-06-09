"""
This is a boilerplate pipeline 'data_collection_pdb'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict, List
import requests
from requests.adapters import HTTPAdapter, Retry
import pandas as pd
import json
import time
from urllib.parse import quote

logger = logging.getLogger(__name__)


def _fetch_pbd_ids() -> List[str]:
    """Fetch the list of PDB ids to be used for data collection"""
    # Fetch the list of PDB ids
    response = requests.get(
        "https://data.rcsb.org/rest/v1/holdings/current/entry_ids", timeout=10
    )
    return response.json()


def _query_pdb(ids: list) -> str:
    ids_string = ", ".join(f'"{id}"' for id in ids)
    return f"""
    {{
        entries(entry_ids: [{ids_string}]) {{
            rcsb_id
            rcsb_accession_info {{
            initial_release_date
            }}
            audit_author {{
            name
            }}
            refine {{
            ls_R_factor_R_free
            }}
            rcsb_entry_info {{
            resolution_combined
            }}
            rcsb_primary_citation {{
            pdbx_database_id_PubMed
            pdbx_database_id_DOI
            }}
        }}
    }}
    """


def _create_structure_search_request(entry_id: str, assembly_id: str = "1") -> str:
    """
    Create a structure search request JSON for the PDB API.

    Args:
        entry_id (str): The PDB entry ID to search for similar structures.
        assembly_id (str): The assembly ID, defaults to "1".

    Returns:
        str: JSON string for the API request.
    """
    search_request = {
        "query": {
            "type": "terminal",
            "service": "structure",
            "parameters": {
                "value": {"entry_id": entry_id, "assembly_id": assembly_id},
                "operator": "strict_shape_match",
                "target_search_space": "assembly",
            },
        },
        "return_type": "assembly",
        "request_options": {
            "paginate": {"start": 0, "rows": 10},  # Top 10 results as requested
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "scoring_strategy": "combined",
        },
    }

    return json.dumps(search_request)


def _extract_pdb_id(identifier: str) -> str:
    """
    Extract the 4-character PDB ID from an identifier like '1LGA-1'.

    Args:
        identifier (str): Full identifier from API response.

    Returns:
        str: 4-character PDB ID.
    """
    return identifier.split("-")[0][:4]


def _query_structure_matches(
    entry_id: str, session: requests.Session, max_retries: int = 3
) -> pd.DataFrame:
    """
    Query the PDB API for structure matches for a single entry.

    Args:
        entry_id (str): The PDB entry ID to search for.
        session (requests.Session): Configured session for API calls.
        max_retries (int): Maximum number of retries for failed requests.

    Returns:
        pd.DataFrame: DataFrame with columns query, target, score.
    """
    search_request = _create_structure_search_request(entry_id)
    url = f"https://search.rcsb.org/rcsbsearch/v2/query?json={quote(search_request)}"

    for attempt in range(max_retries):
        try:
            response = session.get(url, timeout=30)
            response.raise_for_status()

            data = response.json()

            # Extract results
            if "result_set" not in data or not data["result_set"]:
                logger.warning("No results found for entry %s", entry_id)
                return pd.DataFrame(columns=["query", "target", "score"])

            results = []
            for result in data["result_set"]:
                target_id = _extract_pdb_id(result["identifier"])
                score = result["score"]
                results.append({"query": entry_id, "target": target_id, "score": score})

            df = pd.DataFrame(results)
            logger.info("Found %d structure matches for %s", len(df), entry_id)
            return df

        except requests.exceptions.RequestException as e:
            logger.warning(
                "Attempt %d failed for entry %s: %s", attempt + 1, entry_id, str(e)
            )
            if attempt < max_retries - 1:
                time.sleep(2**attempt)  # Exponential backoff
            else:
                logger.error(
                    "Failed to get structure matches for %s after %d attempts",
                    entry_id,
                    max_retries,
                )
                return pd.DataFrame(columns=["query", "target", "score"])

        except (json.JSONDecodeError, KeyError) as e:
            logger.error("Error parsing response for %s: %s", entry_id, str(e))
            return pd.DataFrame(columns=["query", "target", "score"])


def fetch_pbd_details(config: Dict[str, str]) -> pd.DataFrame:
    """Fetch the details of a PDB id"""

    pbd_ids = _fetch_pbd_ids()
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
        query = _query_pdb(ids)
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

    This function queries the PDB structure search API for each entry to find
    similar structures and their match scores. It returns the top 10 matches
    for each query structure.

    Args:
        pdb_entries (pd.DataFrame): DataFrame containing PDB entries with 'rcsb_id' column.
        config (Dict[str, str]): Configuration parameters including max_retries and backoff_factor.

    Returns:
        pd.DataFrame: DataFrame with columns 'query', 'target', 'score' containing
                     structure match information.
    """
    logger.info(
        "Starting structure match score collection for %d PDB entries", len(pdb_entries)
    )

    # Set up session with retry configuration
    session = requests.Session()
    retries = Retry(
        total=config.get("max_retries", 3),
        backoff_factor=config.get("backoff_factor", 0.3),
        status_forcelist=[429, 500, 502, 503, 504],
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))

    all_matches = []

    # Get list of PDB IDs to query
    pdb_ids = pdb_entries["rcsb_id"].unique()

    for idx, entry_id in enumerate(pdb_ids):
        if idx % 100 == 0:
            logger.info("Processing entry %d/%d: %s", idx + 1, len(pdb_ids), entry_id)

        # Query structure matches for this entry
        matches_df = _query_structure_matches(entry_id, session)

        if not matches_df.empty:
            all_matches.append(matches_df)

        # Add small delay to be respectful to the API
        time.sleep(0.1)

    if all_matches:
        result_df = pd.concat(all_matches, ignore_index=True)
        logger.info("Collected %d structure match records", len(result_df))
    else:
        logger.warning("No structure matches found for any entries")
        result_df = pd.DataFrame(columns=["query", "target", "score"])

    return result_df
