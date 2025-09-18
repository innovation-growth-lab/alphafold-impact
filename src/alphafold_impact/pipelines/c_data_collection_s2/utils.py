"""Utility functions for asynchronous citation data collection from Semantic Scholar.

This module provides the core functionality for making rate-limited API calls to
Semantic Scholar and processing the responses. It handles both forward citations
and references with proper error handling and retries.
"""

import logging
import asyncio
from typing import Sequence, Dict, Union, List, Any
import aiohttp
import pandas as pd
from tqdm import tqdm

logger = logging.getLogger(__name__)


class RateLimiter:
    """
    Rate limiter that ensures minimum delay between operations.
    This of course renders the async version of the function synchronous,
    but we could request higher throughput - at which point we would want async, so
    infrastructure-wise we keep the pipeline async.
    """

    def __init__(self, calls_per_second: float = 1.0):
        self.min_interval = 1.0 / calls_per_second
        self.last_call_time = 0.0
        self._lock = asyncio.Lock()

    async def acquire(self):
        """Wait until we can make another call."""
        async with self._lock:
            now = asyncio.get_event_loop().time()
            time_since_last_call = now - self.last_call_time
            if time_since_last_call < self.min_interval:
                wait_time = self.min_interval - time_since_last_call
                await asyncio.sleep(wait_time)
            self.last_call_time = asyncio.get_event_loop().time()


async def fetch_citation_details_async(
    work_id: str,
    base_url: str,
    direction: str,
    fields: Sequence[str],
    api_key: str,
    session: aiohttp.ClientSession,
    rate_limiter: RateLimiter,
    perpage: int = 500,
) -> Sequence[Dict[str, str]]:
    """Fetch citation details asynchronously with rate limiting and retries.

    Args:
        work_id (str): The ID of the work to fetch citations for.
        base_url (str): Base URL for the API endpoint.
        direction (str): Direction of citations ('references' or 'citations').
        fields (Sequence[str]): Fields to include in the response.
        api_key (str): API key for authentication.
        session (aiohttp.ClientSession): HTTP session for making requests.
        rate_limiter (RateLimiter): Rate limiter instance to control request frequency.
        perpage (int, optional): Number of results per page. Defaults to 500.

    Returns:
        Sequence[Dict[str, str]]: List of citation details dictionaries.
    """
    offset = 0
    data_list = []
    max_retries = 5
    base_delay = 0.5  # seconds
    retryable_status_codes = {429, 500, 502, 503, 504}

    while True:
        url = (
            f"{base_url}/{work_id}/{direction}?"
            f"fields={','.join(fields)}&offset={offset}&limit={perpage}"
        )

        for attempt in range(max_retries):
            # wait for rate limit
            await rate_limiter.acquire()

            try:
                async with session.get(url, headers={"X-API-KEY": api_key}) as response:
                    if response.status in retryable_status_codes:
                        retry_delay = base_delay * (2**attempt)
                        logger.warning(
                            "Retryable error %d for %s, attempt %d/%d. Waiting %d seconds...",
                            response.status,
                            work_id,
                            attempt + 1,
                            max_retries,
                            retry_delay,
                        )
                        await asyncio.sleep(retry_delay)
                        continue

                    if response.status == 404:
                        logger.debug("Resource not found for %s (404)", work_id)
                        return []
                    elif response.status >= 400:
                        logger.warning(
                            "Non-retryable error %d for %s",
                            response.status,
                            work_id,
                        )
                        return []

                    response.raise_for_status()
                    data = await response.json()
                    chunk = data.get("data", [])
                    data_list.extend(chunk)

                    if len(chunk) < perpage:
                        return data_list

                    offset += perpage
                    break  # Success, break retry loop

            except aiohttp.ClientError as e:
                if isinstance(e, aiohttp.ClientResponseError):
                    if e.status not in retryable_status_codes:
                        logger.error("Non-retryable error for %s: %s", work_id, str(e))
                        return []

                if attempt == max_retries - 1:  # Last attempt
                    logger.error(
                        "Failed to fetch %s after %d retries: %s",
                        work_id,
                        max_retries,
                        str(e),
                    )
                    return []

                retry_delay = base_delay * (2**attempt)
                logger.warning(
                    "Error fetching %s (attempt %d/%d): %s. Retrying in %d seconds...",
                    work_id,
                    attempt + 1,
                    max_retries,
                    str(e),
                    retry_delay,
                )
                await asyncio.sleep(retry_delay)

    return data_list


async def iterate_citation_detail_points_async(
    oa: str,
    parent_doi: str,
    doi: str,
    parent_pmid: str,
    pmid: str,
    direction: str,
    session: aiohttp.ClientSession,
    pbar: tqdm,
    rate_limiter: RateLimiter,
    **kwargs,
) -> Dict[str, Union[str, Sequence[Dict[str, Any]]]]:
    """Async version of iterate_citation_detail_points with rate limiting.

    Args:
        oa (str): The ID of the work to fetch citations for.
        parent_doi (str): The DOI of the parent work.
        doi (str): The DOI of the child work.
        parent_pmid (str): The PMID of the parent work.
        pmid (str): The PMID of the child work.
        direction (str): Direction of citations ('references' or 'citations').
        session (aiohttp.ClientSession): HTTP session for making requests.
        pbar (tqdm): Progress bar for tracking progress.
        rate_limiter (RateLimiter): Rate limiter instance to control request frequency.
        **kwargs: Additional keyword arguments.

    Returns:
        Dict[str, Union[str, Sequence[Dict[str, Any]]]]: Dictionary containing
            the citation details.

    Raises:
        Exception: If there is an error fetching the citation details.
    """
    logger.debug("Fetching citation details for %s", oa)

    for prefix, id_ in [("DOI:", doi), ("PMID:", pmid)]:
        if not id_:
            logger.debug("No relevant %s id for %s", prefix[:-1], oa)
            continue

        work_id = f"{prefix}{id_}"
        try:
            data = await fetch_citation_details_async(
                work_id=work_id,
                direction=direction,
                session=session,
                rate_limiter=rate_limiter,
                **kwargs,
            )

            if direction == "references":
                for item in data:
                    if (eids := item["citedPaper"]["externalIds"]) is not None:
                        if (
                            eids.get("DOI", "") == parent_doi
                            or eids.get("PubMed", "") == parent_pmid
                        ):
                            logger.debug("Found citation details using %s", work_id)
                            pbar.update(1)
                            return item
                logger.debug("No relevant citation found in %s", work_id)
            else:
                logger.debug("Found citation details using %s", work_id)
                pbar.update(1)
                return data

        except Exception as e:  # pylint: disable=W0718
            logger.error("Failed to fetch citation details for %s: %s", work_id, str(e))

    pbar.update(1)
    return {}


async def get_intent_level_0_async(
    oa_dataset: pd.DataFrame,
    **kwargs,
) -> pd.DataFrame:
    """Async version of get_intent_level_0 - queries references from each child paper.

    Args:
        oa_dataset (pd.DataFrame): The dataset containing the works to process.
        **kwargs: Additional keyword arguments.

    Returns:
        pd.DataFrame: The processed dataset containing the references.
    """
    inputs = (
        oa_dataset[oa_dataset["level"] == 0]
        .apply(
            lambda x: (x["id"], x["parent_doi"], x["doi"], x["parent_pmid"], x["pmid"]),
            axis=1,
        )
        .tolist()
    )

    rate_limiter = RateLimiter(calls_per_second=1.0)
    connector = aiohttp.TCPConnector(limit=5)
    timeout = aiohttp.ClientTimeout(total=None)

    async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
        pbar = tqdm(total=len(inputs), desc="Processing Level 0 (references)")
        tasks = [
            iterate_citation_detail_points_async(
                *input,
                direction="references",
                session=session,
                pbar=pbar,
                rate_limiter=rate_limiter,
                **kwargs,
            )
            for input in inputs
        ]
        level_outputs = await asyncio.gather(*tasks)
        pbar.close()

    # Use the reference processing for level 0
    level_dict = dict(
        list(zip(oa_dataset[oa_dataset["level"] == 0]["doi"], level_outputs))
    )
    processed_references = process_intent_references(level_dict)

    return pd.DataFrame(
        processed_references,
        columns=[
            "parent_doi",
            "parent_pmid",
            "doi",
            "influential",
            "intent",
        ],
    )


async def get_intent_level_n_async(
    oa_dataset: pd.DataFrame,
    level: Union[int, None],
    **kwargs,
) -> pd.DataFrame:
    """Async version of get_intent_level for levels > 0 - queries forward citations.

    Args:
        oa_dataset (pd.DataFrame): The dataset containing the works to process.
        level (Union[int, None]): The level of the intent to process. If None, process all levels.
        **kwargs: Additional keyword arguments.

    Returns:
        pd.DataFrame: The processed dataset containing the forward citations.
    """
    level_data = (
        oa_dataset if level is None else oa_dataset[oa_dataset["level"] == level]
    )
    level_data = level_data.drop_duplicates(subset="parent_id")

    inputs = level_data.apply(
        lambda x: (x["parent_id"], "", x["parent_doi"], "", x["parent_pmid"]),
        axis=1,
    ).tolist()

    rate_limiter = RateLimiter(calls_per_second=1.0)
    connector = aiohttp.TCPConnector(limit=5)
    timeout = aiohttp.ClientTimeout(total=None)

    async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
        pbar = tqdm(
            total=len(inputs),
            desc=(
                f"Processing {'All Levels' if level is None else f'Level {level}'} "
                "(forward citations)"
            ),
        )
        tasks = [
            iterate_citation_detail_points_async(
                *input,
                direction="citations",
                session=session,
                pbar=pbar,
                rate_limiter=rate_limiter,
                **kwargs,
            )
            for input in inputs
        ]
        level_outputs = await asyncio.gather(*tasks)
        pbar.close()

    level_dict = dict(list(zip(level_data["parent_doi"], level_outputs)))
    processed_level_citations = process_intent_citations(level_dict)

    return pd.DataFrame(
        processed_level_citations,
        columns=[
            "parent_doi",
            "pmid",
            "doi",
            "influential",
            "intent",
        ],
    )


def process_intent_references(reference_outputs: Dict[str, Dict[str, Any]]) -> List:
    """
    Process the reference outputs with simplified intents field and generate a list of rows
    containing relevant information.

    Args:
        reference_outputs (Dict[str, Dict[str, Any]]): A dictionary containing the
            reference outputs with simplified intents field.

    Returns:
        List: A list of rows containing the processed information with format:
            [parent_doi, parent_pmid, child_doi, influential, intent]
    """
    rows = []
    for child_doi, output in reference_outputs.items():
        influential = output.get("isInfluential", "")
        intents = output.get("intents", [])  # Get the list of intents directly
        parent_doi = output.get("citedPaper", {}).get("externalIds", {}).get("DOI", "")
        parent_pmid = (
            output.get("citedPaper", {}).get("externalIds", {}).get("PubMed", "")
        )

        # If no intents are found, add a row with empty intent
        if not intents:
            rows.append([parent_doi, parent_pmid, child_doi, influential, ""])
        else:
            # Add a row for each intent in the list
            for intent in intents:
                rows.append([parent_doi, parent_pmid, child_doi, influential, intent])

    return rows


def process_intent_citations(
    citation_outputs: Dict[str, Sequence[Dict[str, Any]]],
) -> List:
    """
    Process citation outputs and extract relevant information.

    Args:
        citation_outputs (Dict[str, Sequence[Dict[str, Any]]]): A dictionary
            containing citation outputs.

    Returns:
        List: A list of rows containing extracted information from the
            citation outputs.
    """
    rows = []
    for parent_doi, outputs in citation_outputs.items():
        for output in outputs:
            intent_groups = output.get("intents", [])
            influential = output.get("isInfluential", "")
            external_ids = output.get("citingPaper", {}).get("externalIds", {})
            if not external_ids:
                continue
            child_doi = external_ids.get("DOI", "")
            child_pmid = external_ids.get("PubMed", "")
            for intent in intent_groups:
                rows.append(
                    [
                        parent_doi,
                        child_pmid,
                        child_doi,
                        influential,
                        intent,
                    ]
                )

    return rows


def create_parent_child_df(oa_dataset: pd.DataFrame) -> pd.DataFrame:
    """
    Create a parent-child dataframe from the oa_dataset.

    Args:
        oa_dataset (pd.DataFrame): The dataset containing the works to process.

    Returns:
        pd.DataFrame: The parent-child dataframe.
    """
    inner_df = oa_dataset.copy()
    # if level is str, convert to int with "seed" as -1
    inner_df["level"] = inner_df["level"].apply(
        lambda x: -1 if x == "seed" else int(x)
    )

    # drop duplicates within level
    inner_df = inner_df.drop_duplicates(subset=["id", "parent_id", "level"])

    # Create a mapping from id to doi and pmid for each level
    inner_df["parent_level"] = inner_df["level"] - 1

    inner_df = pd.merge(
        inner_df,
        inner_df,
        left_on=["parent_id", "parent_level"],
        right_on=["id", "level"],
        how="left",
        suffixes=("", "_parent"),
    )
    inner_df.rename(
        columns={"doi_parent": "parent_doi", "pmid_parent": "parent_pmid"}, inplace=True
    )

    # drop the other _parent columns
    inner_df.drop(
        columns=[col for col in inner_df.columns if "_parent" in col], inplace=True
    )

    # create a dictionary to map parent_id to DOI and PMID
    parent_info = {
        "W3177828909": {"doi": "10.1038/s41586-021-03819-2", "pmid": "34265844"},
        "W3211795435": {"doi": "10.1093/nar/gkab1061", "pmid": "34791371"},
        "W3202105508": {"doi": "10.1101/2021.10.04.463034", "pmid": None},
    }

    # extend parent_info with all entries that have level -1 (seed level)
    if -1 in inner_df["level"].values:
        seed_entries = inner_df[inner_df["level"] == -1]
        for _, row in seed_entries.iterrows():
            parent_info[str(row["id"])] = {
                "doi": row["doi"] if pd.notna(row["doi"]) else None,
                "pmid": row["pmid"] if pd.notna(row["pmid"]) else None,
            }

    # replace the values in the 'parent_doi' and 'parent_pmid' columns when level is 0
    inner_df.loc[inner_df["level"] == 0, "parent_doi"] = inner_df.loc[
        inner_df["level"] == 0, "parent_id"
    ].map(lambda x: parent_info[str(x)]["doi"])
    inner_df.loc[inner_df["level"] == 0, "parent_pmid"] = inner_df.loc[
        inner_df["level"] == 0, "parent_id"
    ].map(lambda x: parent_info[str(x)]["pmid"])

    return inner_df
