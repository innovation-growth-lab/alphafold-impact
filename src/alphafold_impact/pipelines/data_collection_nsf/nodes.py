"""
This module contains the functions used to fetch data from the NSF API.

Functions:
    _create_publication_objects: Create publication objects from a list of paper data.
    fetch_award_data: Fetch award data for a given award ID.
    fetch_nsf_data: Fetch NSF data for a given award year.
"""
import logging
from typing import Dict, List
import requests
from requests.adapters import HTTPAdapter, Retry
from urllib3.exceptions import MaxRetryError
from joblib import Parallel, delayed

logger = logging.getLogger(__name__)


def _create_publication_objects(publication_list: List[str]) -> Dict[str, str]:
    """
    Create publication objects from a list of paper data.

    Args:
        publication_list (List[str]): A list of paper data in the format 
            "journal~year~volume~authors~doi~...".

    Returns:
        Dict[str, str]: A dictionary of publication objects, where the keys 
            are the NSF PAR IDs and the values are dictionaries
        containing the publication information (journal, year, volume, 
            authors, title, doi).
    """
    publication_dict = {}
    for paper in publication_list:
        paper_data = paper.split("~")
        if len(paper_data) != 17:
            logger.info("Paper data is incomplete. Skipping paper.")
            continue
        else:
            journal = paper_data[0]
            year = paper_data[1]
            volume = paper_data[2]
            authors = paper_data[3]
            doi = paper_data[4]
            title = paper_data[9]
            nsf_par_id = paper_data[13]

            publication_dict[nsf_par_id] = {
                "journal": journal,
                "year": year,
                "volume": volume,
                "authors": authors,
                "title": title,
                "doi": doi,
            }

    return publication_dict


def fetch_award_data(
    config: Dict[str, str], fields: List[str], award_id: str
) -> Dict[str, str]:
    """
    Fetches award data for a given award ID.

    Args:
        config (Dict[str, str]): Configuration parameters.
        fields (List[str]): List of fields to include in the response.
        award_id (str): ID of the award to fetch data for.

    Returns:
        Dict[str, str]: Award data as a dictionary.

    """
    base_url = config["base_url"]
    fields = ",".join(fields)
    session = requests.Session()
    retries = Retry(
        total=config["max_retries"], backoff_factor=config["backoff_factor"]
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))

    logger.info("Requesting award data for %s", award_id)
    try:
        response = session.get(
            f"{base_url}?id={award_id}?printFields={fields}",
            timeout=30,
        )

        data = response.json()
        award = data.get("response", {}).get("award", [])
        if award:
            logger.info("Received award data for %s", award_id)
            award_object = award[0]
            award_object["publicationResearch"] = _create_publication_objects(
                award_object["publicationResearch"]
            )
            return award_object
        else:
            logger.info("No award data received for %s", award_id)
            return {}
    except MaxRetryError:
        logger.exception("Failed to fetch award data for %s", award_id)
        return {}


def fetch_nsf_data(config: Dict[str, str], fields: List[str], award_year: str):
    """
    Fetches NSF data for a given award year.

    Args:
        config (Dict[str, str]): Configuration parameters for fetching data.
        fields (List[str]): List of fields to fetch for each award.
        award_year (str): Year of the awards to fetch.

    Returns:
        dict: A dictionary containing the fetched NSF data.
    """
    award_ids = [f"{award_year}{i:05}" for i in range(100000)]
    awards = Parallel(n_jobs=6, verbose=10)(
        delayed(fetch_award_data)(config, fields, award_id) for award_id in award_ids
    )
    return {k: v for d in awards for k, v in d.items() if d}
