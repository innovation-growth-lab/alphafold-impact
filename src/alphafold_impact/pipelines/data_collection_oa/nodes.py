import requests
import math
from typing import List
import logging
import pandas as pd


logger = logging.getLogger(__name__)


def create_concept_year_filter(concept_ids: list, years: list) -> str:
    """
    Creates an API query filter string for the OpenAlex API to retrieve works based on a list of concept IDs and years.

    Args:
        concept_ids: A list of concept IDs (e.g., ['c12345', 'c67890']).
        years: A list of publication years (e.g., [2020, 2021]).

    Returns:
        A formatted API query filter string.

    Example:
        >>> create_openalex_filter(['c12345', 'c67890'], [2020, 2021])
        'publication_year:2020|2021,concepts.id:c12345|c67890'
    """
    year_filter = "|".join(map(str, years))  # e.g., "2020|2021"
    if year_filter:
        year_filter = f"publication_year:{year_filter}"

    concept_filter = "|".join(concept_ids)  # e.g., "c12345|c67890"
    if concept_filter:
        concept_filter = f"concepts.id:{concept_filter}"

    # Combine and return full filter string
    return ",".join(
        filter(None, [year_filter, concept_filter])
    )  # Removes empty filters


def retrieve_oa_works_chunk(
    concept_ids: List[str], years: List[int], works_per_page: int
) -> List[dict]:
    """
    Retrieves OpenAlex works for lists of concept IDs and publication years.

    Args:
        concept_ids: OpenAlex concept IDs.
        years: Publication years.
        works_per_page: The number of results to retrieve per page.

    Returns:
        List[dict]: A list of works for the given concept IDs and years.
    """
    base_url = "https://api.openalex.org/works"
    works = []
    next_cursor = "*"
    current_page = 1

    while True:
        params = {
            "filter": create_concept_year_filter(concept_ids, years),
            "per_page": works_per_page,
            "cursor": next_cursor,
        }
        response = requests.get(base_url, params=params)
        if response.status_code == 200:
            data = response.json()
            current_works = data.get("results", [])
            if current_page == 1:
                # Calculate the total number of pages from the first response
                total_count = data.get("meta", {}).get("count", 0)
                total_pages = math.ceil(total_count / works_per_page)

            if not current_works:
                break

            logger.info(f"Processing page {current_page}/{total_pages}")
            works.extend(current_works)
            current_page += 1

            # Retrieve the next cursor from the response
            next_cursor = data.get("meta", {}).get("next_cursor")
            if not next_cursor:
                break
        else:
            logger.error(
                f"Error fetching data for page {current_page}: {response.status_code}"
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
    The same works could be collected more than once if they are related to a concept
    in more than one chunk. Therefore, the works are deduplicated.

    Args:
        concept_ids (List[str]): A list of OpenAlex concept IDs.
        publication_years (List[int]): A list of publication years.
        chunk_size (int): The number of concept IDs to process in each API call (max/default is 40).
        per_page (int): The number of results to retrieve per page (max/default is 200).

    Returns:
        List[dict]: A list of all works that match the criteria.
    """
    all_works = []
    concept_id_chunks = chunk_list(concept_ids, chunk_size)
    total_chunks = len(concept_id_chunks)

    for chunk_index, concept_chunk in enumerate(concept_id_chunks, start=1):
        logger.info(f"Processing concept chunk {chunk_index}/{total_chunks}")
        works = retrieve_oa_works_chunk(concept_chunk, publication_years, per_page)
        all_works.extend(works)
    n_all_works = len(all_works)
    logger.info(f"Retrieved {n_all_works} works.")
    deduped_works = pd.DataFrame(all_works).drop_duplicates(subset=["id"])
    n_deduped_works = len(deduped_works)
    n_dropped = n_all_works - n_deduped_works
    logger.info(f"Dropped {n_dropped} duplicate works.")
    logger.info(f"Saving {n_deduped_works} works.")
    return deduped_works
