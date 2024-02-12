"""Pipeline nodes for the data_collection_s2 pipeline.

Functions:
    load_alphafold_citation_ids: Loads the citation IDs for the AlphaFold paper.
    fetch_citation_details: Fetches citation details for a given work ID.
    get_citation_details: Retrieves citation details for a given list of work IDs 
        and filters them based on a specified AlphaFold DOI.
"""

import logging
from typing import Sequence, Dict, Tuple
import re
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from kedro.io import AbstractDataset

logger = logging.getLogger(__name__)


def load_alphafold_citation_ids(
    input_loaders: AbstractDataset, work_id: str
) -> Sequence[Tuple[str]]:
    """Loads the citation IDs for the AlphaFold paper. Requires OA pipeline."""

    logger.info("Loading citation IDs for Alphafold")
    citations = input_loaders[work_id]()

    # iterate over four possible ids (oa, doi, mag, pmid)
    ids = []
    for citation in citations:
        id_dict = citation.get("ids", {})
        openalex = id_dict.get("openalex", "").replace("https://openalex.org/", "")
        doi = re.sub(r".*?(10\..*)", "\\1", id_dict.get("doi", ""))
        mag = id_dict.get("mag", "")
        pmid = (
            id_dict.get("pmid", "").split("/")[-1]
            if "/" in id_dict.get("pmid", "")
            else id_dict.get("pmid", "")
        )
        ids.append((openalex, doi, mag, pmid))

    return ids


def fetch_citation_details(
    work_id: str,
    base_url: str,
    direction: str,
    fields: Sequence[str],
    perpage: int = 500,
) -> Sequence[Dict[str, str]]:
    """
    Fetches citation details for a given work ID.

    Args:
        work_id (str): The work ID to fetch citation details for.
        base_url (str): The base URL for the API.
        direction (str): The direction of the citations.
        fields (Sequence[str]): The fields to fetch.
        perpage (int, optional): The number of citations to fetch per page.
            Defaults to 500.

    Returns:
        Sequence[Dict[str, str]]: A list of citation details.
    """
    offset = 0
    url = (
        f"{base_url}/{work_id}/{direction}?"
        f"fields={','.join(fields)}&offset={offset}&limit={perpage}"
    )

    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=0.3,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))

    data_list = []

    while True:
        response = session.get(url)
        response.raise_for_status()
        data = response.json().get("data", [])
        data_list.extend(data)

        if len(data) < perpage:
            break

        offset += perpage
        url = (
            f"{base_url}/{work_id}/{direction}?"
            f"fields={','.join(fields)}&offset={offset}&limit={perpage}"
        )

    return data_list


def get_alphafold_citation_details(
    work_ids: Sequence[Tuple[str, str, str, str]], af_doi: str, **kwargs
) -> Dict[str, Dict[str, str]]:
    """Retrieves citation details for a given list of work IDs and filters
    them based on the specified AlphaFold DOI.

    Args:
        work_ids (Sequence[Tuple[str, str, str, str]]): A list of work IDs.
        af_doi (str): The AlphaFold DOI to filter the citation details.
        **kwargs: Additional keyword arguments to be passed to the
            fetch_citation_details function.

    Returns:
        Dict[str, Dict[str, str]]: A dictionary containing the filtered
            citation details, where the keys are the work IDs and the values
            are the corresponding citation details.
    """
    citation_details = {}
    for work_id_tuple in work_ids:
        oa, doi, mag, pmid = work_id_tuple
        logger.info("Fetching citation details for %s", oa)
        for prefix, id_ in [("DOI:", doi), ("PMID:", pmid), ("MAG:", mag)]:
            if not id_:
                logger.warning("No relevant %s id", prefix[:-1])
                continue
            work_id = f"{prefix}{id_}"
            try:
                data = fetch_citation_details(work_id, **kwargs)
                for item in data:
                    if (eids := item["citedPaper"]["externalIds"]) is not None:
                        if eids.get("DOI", "") == af_doi:
                            logger.info("Found citation details using %s", work_id)
                            citation_details[oa] = item
                            break
                else:
                    logger.warning("No relevant citation found in %s", work_id)
                    continue
                break
            except Exception:  # pylint: disable=broad-except
                logger.warning("Failed to fetch citation details for %s", work_id)

    # explode multiple citations
    rows = []
    for child, details in citation_details.items():
        intents = details.get("intents", [])
        contexts = details.get("contexts", [])
        for intent, context in zip(intents, contexts):
            rows.append(["W3177828909", child, intent, context])

    # create dataframe
    citation_details = pd.DataFrame(
        rows, columns=["parent", "child", "intent", "context"]
    )

    return citation_details
