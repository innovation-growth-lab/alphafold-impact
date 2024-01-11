"""
This module contains functions for fetching and preprocessing data from the OpenAlex API.

Functions:
    collect_citation_papers: Collects all papers cited by specific work IDs.
    load_work_ids: Loads the file corresponding to a particular work_id in a PartitionedDataSet,
        extracts all ids, and returns these as a list.
    retrieve_oa_works_for_concepts_and_years: Retrieves OpenAlex works for the specified
        concept IDs and publication years.

Internal functions:
    _revert_abstract_index: Reverts the abstract inverted index to the original text.
    _parse_results: Parses OpenAlex API response to retain basic variables.
    _citation_works_generator: Creates a generator that yields a list of works from the OpenAlex API based on a
        given work ID.
    _create_concept_year_filter: Creates an API query filter string for the OpenAlex API.
    _retrieve_oa_works_chunk: Retrieves a chunk of OpenAlex works for a chunk of concept IDs.
    _chunk_list: Divides the input list into chunks of specified size.
"""
import logging
from typing import Iterator, List, Dict, Sequence, Union
from requests.adapters import HTTPAdapter, Retry
import requests
from kedro.io import AbstractVersionedDataset
import pandas as pd


logger = logging.getLogger(__name__)


def _revert_abstract_index(abstract_inverted_index: Dict[str, Sequence[int]]) -> str:
    """Reverts the abstract inverted index to the original text.

    Args:
        abstract_inverted_index (Dict[str, Sequence[int]]): The abstract inverted index.

    Returns:
        str: The original text.
    """
    try:
        length_of_text = (
            max(
                [
                    index
                    for sublist in abstract_inverted_index.values()
                    for index in sublist
                ]
            )
            + 1
        )
        recreated_text = [""] * length_of_text

        for word, indices in abstract_inverted_index.items():
            for index in indices:
                recreated_text[index] = word

        return " ".join(recreated_text)
    except (AttributeError, ValueError):
        return ""


def _parse_results(response: List[Dict]) -> Dict[str, List[str]]:
    """Parses OpenAlex API response to retain:
        - id
        - display_name
        - title
        - publication_date
        - abstract_inverted_index
        - authorships
        - cited_by_count
        - concepts
        - keywords
        - grants

    Args:
        response (List[Dict]): The response from the OpenAlex API.

    Returns:
        Dict[str, List[str]]: A dictionary containing the parsed information.
    """
    return [
        {
            "id": paper["id"],
            "display_name": paper["display_name"],
            "title": paper["title"],
            "publication_date": paper["publication_date"],
            "abstract": _revert_abstract_index(paper["abstract_inverted_index"]),
            "authorships": paper["authorships"],
            "cited_by_count": paper["cited_by_count"],
            "concepts": paper["concepts"],
            "keywords": paper["keywords"],
            "grants": paper["grants"],
        }
        for paper in response
    ]


def _citation_works_generator(
    mailto: str, perpage: str, work_id: str, direction: str
) -> Iterator[list]:
    """Creates a generator that yields a list of works from the OpenAlex API based on a
    given work ID.

    Args:
        mailto (str): The email address to use for the API.
        perpage (str): The number of results to return per page.
        work_id (str): A single work ID to filter by 'cites'.
        direction (str): The direction of the citation. Either 'cites' or 'cited_by'.

    Yields:
        Iterator[list]: A generator that yields a list of works from the OpenAlex API
        based on a given work ID.
    """

    cursor_url = (
        f"https://api.openalex.org/works?filter={direction}:{work_id}"
        f"&mailto={mailto}&per-page={perpage}&cursor={{}}"
    )

    cursor = "*"
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=0.3)
    session.mount("https://", HTTPAdapter(max_retries=retries))

    # make a call to estimate total number of results
    response = session.get(cursor_url.format(cursor), timeout=20)
    data = response.json()
    total_results = data["meta"]["count"]
    num_calls = total_results // int(perpage) + 1
    logger.info(
        "Total results for %s: %s, in %s calls", work_id, total_results, num_calls
    )

    while cursor:
        response = session.get(cursor_url.format(cursor), timeout=20)
        data = response.json()
        results = data.get("results")
        cursor = data["meta"].get("next_cursor", False)
        yield results


def collect_citation_papers(
    mailto: str, perpage: str, work_ids: Union[str, List[str]], direction: str
) -> dict:
    """Collects all papers cited by specific work IDs.

    Args:
        mailto (str): The email address to use for the API.
        work_ids (List[str]): A list of work IDs to filter by either 'cites' or 'cited_by'.
        direction (str): The direction of the citation. Either 'cites' or 'cited_by'.
        perpage (str): The number of results to return per page.

    Returns:
        dict: A dictionary containing the work_id as key and the list of papers as value.
    """

    assert direction in ["cites", "cited_by"], (
        f"Invalid direction: {direction}. " f"Must be either 'cites' or 'cited_by'."
    )

    if isinstance(work_ids, str):
        work_ids = [work_ids]

    all_papers = {}
    for work_id in work_ids:
        papers_for_id = []
        for page, papers in enumerate(
            _citation_works_generator(
                mailto=mailto, perpage=perpage, work_id=work_id, direction=direction
            )
        ):
            papers_for_id.extend(_parse_results(papers))
            logger.info(
                "Fetching page %s of %s. Total papers collected: %s",
                page,
                work_id,
                len(papers_for_id),
            )
        all_papers[work_id] = papers_for_id

    return all_papers


def load_work_ids(
    work_id: str, dataset: Sequence[AbstractVersionedDataset]
) -> List[str]:
    """
    Loads the file corresponding to a particular work_id in a PartitionedDataSet,
    extracts all ids, and returns these as a list.

    Args:
        work_id (str): The work_id to load the file for.
        dataset (PartitionedDataSet): The PartitionedDataSet containing the files.

    Returns:
        List[str]: A list of all ids extracted from the file.
    """
    data = dataset[work_id]()
    ids = [paper["id"].replace("https://openalex.org/", "") for paper in data]
    logger.info("Loaded %s ids for %s", len(ids), work_id)
    return ids


def _create_concept_year_filter(concept_ids: List[str], years: List[int]) -> str:
    """
    Creates an API query filter string for the OpenAlex API to retrieve works
    based on lists of concept IDs and years.

    Args:
        concept_ids (List[str]): A list of concept IDs (e.g., ['c12345', 'c67890']).
        years (List[int]): A list of publication years (e.g., [2020, 2021]).

    Returns:
        str: A formatted API query filter string.
    """
    year_filter = f"publication_year:{'|'.join(map(str, years))}" if years else ""
    concept_filter = f"concepts.id:{'|'.join(concept_ids)}" if concept_ids else ""
    return ",".join(filter(None, [year_filter, concept_filter]))


def _retrieve_oa_works_chunk(
    concept_ids: List[str], years: List[int], works_per_page: int
) -> List[dict]:
    """
    Retrieves a chunk of OpenAlex works for a chunk of concept IDs
    and publication years.

    Args:
        concept_ids (List[str]): A list of OpenAlex concept IDs.
        years (List[int]): A list of publication years.
        works_per_page (int): The number of results to retrieve per API request.

    Returns:
        List[dict]: List containing works for the given concept IDs and years.
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
            "filter": _create_concept_year_filter(concept_ids, years),
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
                logger.error("Error fetching data: %s", response.status_code)
                break
        except requests.exceptions.RequestException as e:
            logger.error("Request failed: %s", e)
            break

    return works


def _chunk_list(input_list: List[str], chunk_size: int) -> List[List[str]]:
    """
    Divides the input list into chunks of specified size.

    Args:
        input_list (List[str]): The input list to be divided into chunks.
        chunk_size (int): The size of each chunk.

    Returns:
        List[List[str]]: A list of chunks.
    """
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
        concept_ids (List[str]): A list of OpenAlex concept IDs.
        publication_years (List[int]): A list of publication years.
        chunk_size (int, optional): The number of concept IDs to process in each API call (default is 40).
        per_page (int, optional): The number of results to retrieve per API call (default is 200).

    Returns:
        pd.DataFrame: A pandas DataFrame containing all unique retrieved works.
    """
    all_works = []
    concept_id_chunks = _chunk_list(concept_ids, chunk_size)

    for index, concept_chunk in enumerate(concept_id_chunks, start=1):
        logger.info("Processing chunk %s/%s", index, len(concept_id_chunks))
        chunk_works = _retrieve_oa_works_chunk(
            concept_chunk, publication_years, per_page
        )
        all_works.extend(chunk_works)

    all_works_df = pd.DataFrame(all_works).drop_duplicates(subset=["id"])
    logger.info("Retrieved %s works.", len(all_works_df))
    return all_works_df
