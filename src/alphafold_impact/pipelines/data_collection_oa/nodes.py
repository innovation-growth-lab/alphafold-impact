import requests
import math
from typing import List
from tqdm import tqdm
import pandas as pd


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
                page_progress = tqdm(
                    total=total_pages, desc="Pages", leave=False
                )  # Initialize progress bar with total
                first_request = False

            if not current_works:  # Break if there are no more works
                break
            works.extend(current_works)

            # Update the page progress bar
            page_progress.update(1)

            # Retrieve the next cursor from the response
            next_cursor = data.get("meta", {}).get("next_cursor")
            if not next_cursor:  # Break if there's no next cursor
                break
        else:
            print(
                f"Error fetching data for concept {concept_id} and year {year}: {response.status_code}"
            )
            break

    if not first_request:
        page_progress.close()

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

    for year in tqdm(publication_years, desc="Years"):
        for concept_chunk in tqdm(concept_chunks, desc="Concept Chunks", leave=False):
            for concept_id in tqdm(concept_chunk, desc="Concepts", leave=False):
                works = retrieve_oa_works_for_single_concept_and_year(
                    concept_id, year, per_page
                )
                all_works.extend(works)
    n_works = len(works)
    print(f"Retrieved {n_works} works")
    return works


def save_oa_works_to_s3(works: List[dict]) -> JSONDataset:
    """Transforms OpenAlex works list into json object and saves to S3"""
    json_dataset = JSONDataset(
        filepath="s3://alphafold-impact/data/01_raw/openalex/works/works_for_concepts.json"
    )
    json_dataset.save(works)
    return json_dataset
