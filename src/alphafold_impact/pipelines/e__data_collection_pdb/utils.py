"""Utility functions for PDB data collection pipeline.

This module contains helper functions for:
- Fetching PDB IDs and details
- Creating structure search requests
- Parsing API responses
- Asynchronous structure matching with rate limiting and retries
"""

import logging
import asyncio
import json
from typing import List
from urllib.parse import quote

import aiohttp
import pandas as pd
import requests
from tqdm import tqdm
from ..c_data_collection_s2.utils import RateLimiter  # pylint: disable=E0402

logger = logging.getLogger(__name__)


def fetch_pdb_ids() -> List[str]:
    """Fetch the list of PDB ids to be used for data collection"""
    # Fetch the list of PDB ids
    response = requests.get(
        "https://data.rcsb.org/rest/v1/holdings/current/entry_ids", timeout=10
    )
    return response.json()


def create_pdb_query(ids: list) -> str:
    """Create GraphQL query for PDB entries.

    Args:
        ids (list): List of PDB IDs to query.

    Returns:
        str: GraphQL query string.
    """
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


def create_structure_search_request(entry_id: str, assembly_id: str = "1") -> str:
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


def extract_pdb_id(identifier: str) -> str:
    """
    Extract the 4-character PDB ID from an identifier like '1LGA-1'.

    Args:
        identifier (str): Full identifier from API response.

    Returns:
        str: 4-character PDB ID.
    """
    return identifier.split("-")[0][:4]


async def query_structure_matches_async(
    entry_id: str,
    session: aiohttp.ClientSession,
    rate_limiter: RateLimiter,
    pbar: tqdm,
    max_retries: int = 5,
) -> pd.DataFrame:
    """
    Query the PDB API for structure matches for a single entry asynchronously.

    Args:
        entry_id (str): The PDB entry ID to search for.
        session (aiohttp.ClientSession): HTTP session for making requests.
        rate_limiter (RateLimiter): Rate limiter instance to control request frequency.
        pbar (tqdm): Progress bar for tracking progress.
        max_retries (int): Maximum number of retries for failed requests.

    Returns:
        pd.DataFrame: DataFrame with columns query, target, score.
    """
    search_request = create_structure_search_request(entry_id)
    url = f"https://search.rcsb.org/rcsbsearch/v2/query?json={quote(search_request)}"

    base_delay = 0.5  # seconds
    retryable_status_codes = {429, 500, 502, 503, 504}

    for attempt in range(max_retries):
        # Wait for rate limit
        await rate_limiter.acquire()

        try:
            async with session.get(url) as response:
                if response.status in retryable_status_codes:
                    retry_delay = base_delay * (2**attempt)
                    logger.warning(
                        "Retryable error %d for %s, attempt %d/%d. Waiting %.1f seconds...",
                        response.status,
                        entry_id,
                        attempt + 1,
                        max_retries,
                        retry_delay,
                    )
                    await asyncio.sleep(retry_delay)
                    continue

                if response.status == 404:
                    logger.debug("No results found for entry %s (404)", entry_id)
                    pbar.update(1)
                    return pd.DataFrame(columns=["query", "target", "score"])
                elif response.status >= 400:
                    logger.warning(
                        "Non-retryable error %d for %s",
                        response.status,
                        entry_id,
                    )
                    pbar.update(1)
                    return pd.DataFrame(columns=["query", "target", "score"])

                response.raise_for_status()
                data = await response.json()

                # Extract results
                if "result_set" not in data or not data["result_set"]:
                    logger.debug("No results found for entry %s", entry_id)
                    pbar.update(1)
                    return pd.DataFrame(columns=["query", "target", "score"])

                results = []
                for result in data["result_set"]:
                    target_id = extract_pdb_id(result["identifier"])
                    score = result["score"]
                    results.append(
                        {"query": entry_id, "target": target_id, "score": score}
                    )

                df = pd.DataFrame(results)
                logger.debug("Found %d structure matches for %s", len(df), entry_id)
                pbar.update(1)
                return df

        except aiohttp.ClientError as e:
            if isinstance(e, aiohttp.ClientResponseError):
                if e.status not in retryable_status_codes:
                    logger.error("Non-retryable error for %s: %s", entry_id, str(e))
                    pbar.update(1)
                    return pd.DataFrame(columns=["query", "target", "score"])

            if attempt == max_retries - 1:  # Last attempt
                logger.error(
                    "Failed to fetch structure matches for %s after %d retries: %s",
                    entry_id,
                    max_retries,
                    str(e),
                )
                pbar.update(1)
                return pd.DataFrame(columns=["query", "target", "score"])

            retry_delay = base_delay * (2**attempt)
            logger.warning(
                "Error fetching %s (attempt %d/%d): %s. Retrying in %.1f seconds...",
                entry_id,
                attempt + 1,
                max_retries,
                str(e),
                retry_delay,
            )
            await asyncio.sleep(retry_delay)

        except (json.JSONDecodeError, KeyError) as e:
            logger.error("Error parsing response for %s: %s", entry_id, str(e))
            pbar.update(1)
            return pd.DataFrame(columns=["query", "target", "score"])

    # If we reach here, all retries failed
    pbar.update(1)
    return pd.DataFrame(columns=["query", "target", "score"])


async def collect_structure_matches_async(
    pdb_ids: List[str],
    calls_per_second: float = 2.0,
    max_concurrent: int = 10,
    max_retries: int = 5,
) -> pd.DataFrame:
    """
    Collect structure match scores for multiple PDB entries asynchronously.

    Args:
        pdb_ids (List[str]): List of PDB entry IDs to query.
        calls_per_second (float): Rate limit for API calls.
        max_concurrent (int): Maximum number of concurrent requests.
        max_retries (int): Maximum number of retries for failed requests.

    Returns:
        pd.DataFrame: DataFrame with columns query, target, score.
    """
    rate_limiter = RateLimiter(calls_per_second=calls_per_second)
    connector = aiohttp.TCPConnector(limit=max_concurrent)
    timeout = aiohttp.ClientTimeout(total=None)

    async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
        pbar = tqdm(
            total=len(pdb_ids),
            desc="Collecting PDB structure matches",
        )

        tasks = [
            query_structure_matches_async(
                entry_id=pdb_id,
                session=session,
                rate_limiter=rate_limiter,
                pbar=pbar,
                max_retries=max_retries,
            )
            for pdb_id in pdb_ids
        ]

        results = await asyncio.gather(*tasks, return_exceptions=True)
        pbar.close()

    # Filter out exceptions and combine results
    valid_results = []
    for i, result in enumerate(results):
        if isinstance(result, Exception):
            logger.error("Exception for PDB ID %s: %s", pdb_ids[i], str(result))
        elif isinstance(result, pd.DataFrame) and not result.empty:
            valid_results.append(result)

    if valid_results:
        combined_df = pd.concat(valid_results, ignore_index=True)
        logger.info("Collected %d structure match records", len(combined_df))
        return combined_df
    else:
        logger.warning("No structure matches found for any entries")
        return pd.DataFrame(columns=["query", "target", "score"])
