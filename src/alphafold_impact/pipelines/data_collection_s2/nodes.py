"""Pipeline nodes for the data_collection_s2 pipeline.

Functions:
    load_alphafold_citation_ids: Loads the citation IDs for the AlphaFold paper.
    fetch_citation_details: Fetches citation details for a given work ID.
    get_citation_details: Retrieves citation details for a given list of work IDs 
        and filters them based on a specified AlphaFold DOI.
"""

import logging
from typing import Sequence, Dict
import re
import requests
from requests.adapters import HTTPAdapter, Retry
from kedro.io import AbstractDataset

logger = logging.getLogger(__name__)


def load_alphafold_citation_ids(
    input_loaders: AbstractDataset, work_id: str
) -> Sequence[str]:
    """Loads the citation IDs for the AlphaFold paper. Requires OA pipeline."""

    logger.info("Loading citation IDs for Alphafold")
    citations = input_loaders[work_id]()

    dois = [
        re.search(r"10\..*", item["doi"]).group() if item["doi"] else None
        for item in citations
    ]

    # get rid of None values
    dois = [item for item in dois if item]

    return dois


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


def get_citation_details(
    work_ids: Sequence[str], af_doi: str, **kwargs
) -> Dict[str, Dict[str, str]]:
    """Retrieves citation details for a given list of work IDs and filters
    them based on a specified AlphaFold DOI.

    Args:
        work_ids (Sequence[str]): A list of work IDs.
        af_doi (str): The AlphaFold DOI to filter the citation details.
        **kwargs: Additional keyword arguments to be passed to the
            fetch_citation_details function.

    Returns:
        Dict[str, Dict[str, str]]: A dictionary containing the filtered
            citation details, where the keys are the work IDs and the values
            are the corresponding citation details.
    """
    citation_details = {}
    for work_id in work_ids:
        logger.info("Fetching citation details for %s", work_id)
        work_id = f"DOI:{work_id}"
        try:
            data = fetch_citation_details(work_id, **kwargs)
            for item in data:
                if item["citedPaper"]["externalIds"].get("DOI", "") == af_doi:
                    citation_details[work_id] = item
                    break
        except Exception: # pylint: disable=broad-except
            logger.warning("Failed to fetch citation details for %s", work_id)
            continue
    return citation_details

# Some are failing, the result of DOI not being quite unique (the domain does redirects).
# We may be able to fetch these using PubMed IDs, but the baseline OA pipeline did not fetch these.
