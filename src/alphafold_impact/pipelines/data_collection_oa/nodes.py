import requests
import math
from typing import List
import logging
import pandas as pd


logger = logging.getLogger(__name__)


def retrieve_oa_works_for_single_concept_and_year(
    concept_id: str, year: int, per_page: int
) -> List[dict]:
    """
    Retrieves OpenAlex works for a single concept ID and publication year.

    Args:
        concept_id (str): The OpenAlex concept ID.
        year (int): The publication year.
        per_page (int): The number of results to retrieve per page.

    Returns:
        List[dict]: A list of works for the given concept ID and year.
    """
    base_url = "https://api.openalex.org/works"
    works = []
    next_cursor = "*"
    first_request = True

    while True:
        params = {
            "filter": f"concepts.id:{concept_id},publication_year:{year}",
            "per_page": per_page,
            "cursor": next_cursor,
        }
        response = requests.get(base_url, params=params)
        if response.status_code == 200:
            data = response.json()
            current_works = data.get("results", [])
            if first_request:
                # Calculate the total number of pages from the first response
                total_count = data.get("meta", {}).get("count", 0)
                total_pages = math.ceil(total_count / per_page)
                current_page = 1
                first_request = False

            if not current_works:
                break

            logger.info(
                f"Processing concept {concept_id.upper()} page {current_page}/{total_pages}"
            )
            works.extend(current_works)
            current_page += 1

            # Retrieve the next cursor from the response
            next_cursor = data.get("meta", {}).get("next_cursor")
            if not next_cursor:
                break
        else:
            logger.error(
                f"Error fetching data for concept {concept_id} and year {year}: {response.status_code}"
            )
            break

    return works


def retrieve_oa_works_for_single_concept_and_year(
    concept_id: str, year: int, per_page: int
) -> List[dict]:
    """
    Retrieves OpenAlex works for a single concept ID and publication year.

    Args:
        concept_id (str): The OpenAlex concept ID.
        year (int): The publication year.
        per_page (int): The number of results to retrieve per page.

    Returns:
        List[dict]: A list of works for the given concept ID and year.
    """
    base_url = "https://api.openalex.org/works"
    works = []
    next_cursor = "*"
    current_page = 1

    while True:
        params = {
            "filter": f"concepts.id:{concept_id},publication_year:{year}",
            "per_page": per_page,
            "cursor": next_cursor,
        }
        response = requests.get(base_url, params=params)
        if response.status_code == 200:
            data = response.json()
            current_works = data.get("results", [])
            if current_page == 1:
                # Calculate the total number of pages from the first response
                total_count = data.get("meta", {}).get("count", 0)
                total_pages = math.ceil(total_count / per_page)

            if not current_works:
                break  # Break if there are no more works

            logger.info(
                f"Processing page {current_page}/{total_pages} for concept {concept_id.upper()}"
            )
            works.extend(current_works)
            current_page += 1

            # Retrieve the next cursor from the response
            next_cursor = data.get("meta", {}).get("next_cursor")
            if not next_cursor:
                break  # Break if there's no next cursor
        else:
            logger.error(
                f"Error fetching data for concept {concept_id} and year {year}: {response.status_code}"
            )
            break

    return works


def chunk_list(input_list: List[str], chunk_size: int) -> List[List[str]]:
    """Divides the input list into chunks of specified size."""
    return [
        input_list[i : i + chunk_size] for i in range(0, len(input_list), chunk_size)
    ]


def retrieve_oa_works_for_concepts_and_years(
    concept_ids: List[str],
    publication_years: List[int],
    chunk_size: int = 40,
    per_page: int = 200,
) -> List[dict]:
    """
    Retrieves OpenAlex works for the concept IDs and publication years provided.

    Args:
        concept_ids (List[str]): A list of OpenAlex concept IDs.
        publication_years (List[int]): A list of publication years.
        chunk_size (int): The number of concept IDs to process in each API call (max/default is 40).
        per_page (int): The number of results to retrieve per page (max/default is 200).

    Returns:
        List[dict]: A list of all works that match the criteria.
    """
    all_works = []
    concept_chunks = chunk_list(concept_ids, chunk_size)
    total_chunks = len(concept_chunks)

    for year in publication_years:
        logger.info(f"Processing year: {year}")
        for chunk_index, concept_chunk in enumerate(concept_chunks, start=1):
            logger.info(f"Processing concept chunk {chunk_index}/{total_chunks}")
            for concept_id in concept_chunk:
                works = retrieve_oa_works_for_single_concept_and_year(
                    concept_id, year, per_page
                )
                all_works.extend(works)
    n_all_works = len(all_works)
    logger.info(f"Retrieved {n_all_works} works")
    return pd.DataFrame(all_works)
