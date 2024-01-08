import logging
from typing import List
import requests
from requests.adapters import HTTPAdapter, Retry
import pandas as pd


logger = logging.getLogger(__name__)


def create_concept_year_filter(concept_ids: List[str], years: List[int]) -> str:
    """
    Creates an API query filter string for the OpenAlex API to retrieve works
    based on lists of concept IDs and years.

    Args:
        concept_ids: A list of concept IDs (e.g., ['c12345', 'c67890']).
        years: A list of publication years (e.g., [2020, 2021]).

    Returns:
        A formatted API query filter string.
    """
    year_filter = f"publication_year:{'|'.join(map(str, years))}" if years else ""
    concept_filter = f"concepts.id:{'|'.join(concept_ids)}" if concept_ids else ""
    return ",".join(filter(None, [year_filter, concept_filter]))


def retrieve_oa_works_chunk(
    concept_ids: List[str], years: List[int], works_per_page: int
) -> List[dict]:
    """
    Retrieves a chunk of OpenAlex works for a chunk of concept IDs
    and publication years.

    Args:
        concept_ids: A list of OpenAlex concept IDs.
        years: A list of publication years.
        works_per_page: The number of results to retrieve per API request.

    Returns:
        List containing works for the given concept IDs and years.
    """
    base_url = "https://api.openalex.org/works"
    next_cursor = "*"
    works = []

    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=0.3,
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))

    while True:
        params = {
            "filter": create_concept_year_filter(concept_ids, years),
            "per_page": works_per_page,
            "cursor": next_cursor,
        }
        try:
            response = session.get(base_url, params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()
                current_works = data.get("results", [])
                works.extend(current_works)
                next_cursor = data.get("meta", {}).get("next_cursor")
                if not (current_works and next_cursor):
                    break
            else:
                logger.error(f"Error fetching data: {response.status_code}")
                break
        except requests.exceptions.RequestException as e:
            logger.error(f"Request failed: {e}")
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
) -> pd.DataFrame:
    """
    Retrieves OpenAlex works for the specified concept IDs and publication years.

    This function processes the concept IDs in chunks to adhere to API limitations and
    compiles the results into a single deduplicated DataFrame.

    Args:
        concept_ids: A list of OpenAlex concept IDs.
        publication_years: A list of publication years.
        chunk_size: The number of concept IDs to process in each API call (default is 40).
        per_page: The number of results to retrieve per API call (default is 200).

    Returns:
        A pandas DataFrame containing all unique retrieved works.
    """
    all_works = []
    concept_id_chunks = chunk_list(concept_ids, chunk_size)

    for index, concept_chunk in enumerate(concept_id_chunks, start=1):
        logger.info(f"Processing chunk {index}/{len(concept_id_chunks)}")
        chunk_works = retrieve_oa_works_chunk(
            concept_chunk, publication_years, per_page
        )
        all_works.extend(chunk_works)

    all_works_df = pd.DataFrame(all_works).drop_duplicates(subset=["id"])
    logger.info(f"Retrieved {len(all_works_df)} works.")
    return all_works_df
