"""
Utilities for data collection from OpenAlex.
    _revert_abstract_index: Revert the abstract inverted index to the original text.
    _parse_results: Parse OpenAlex API response
    _works_generator: Create a generator that yields a list of works from the OpenAlex API 
        based on a given work ID.
    preprocess_oa_ids: Preprocess oa_ids to ensure they are in the correct format.
    _chunk_oa_ids: Yield successive chunk_size-sized chunks from ids.
    fetch_papers_for_id: Fetch all papers cited by a specific work ID.
    fetch_papers_parallel: Fetch papers in parallel.
    fetch_papers_eager: Fetch papers eagerly.
    fetch_papers_lazy: Fetch papers lazily.
    chunk_list: Divides the input list into chunks of specified size.
"""

import logging
from typing import Iterator, List, Dict, Sequence, Union, Callable, Generator
import time
from requests.adapters import HTTPAdapter, Retry
import requests
from joblib import Parallel, delayed

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
        id, doi, display_name, title, publication_date, abstract, authorships,
            cited_by_count, concepts, keywords, grants, referenced_works

    Args:
        response (List[Dict]): The response from the OpenAlex API.

    Returns:
        Dict[str, List[str]]: A dictionary containing the parsed information.
    """
    return [
        {
            "id": paper.get("id", ""),
            "doi": paper.get("doi", ""),
            "display_name": paper.get("display_name", ""),
            "title": paper.get("title", ""),
            "publication_date": paper.get("publication_date", ""),
            "abstract": _revert_abstract_index(
                paper.get("abstract_inverted_index", {})
            ),
            "authorships": paper.get("authorships", []),
            "cited_by_count": paper.get("cited_by_count", ""),
            "concepts": paper.get("concepts", []),
            "keywords": paper.get("keywords", []),
            "mesh_terms": paper.get("mesh", []),
            "topics": paper.get("topics", []),
            "grants": paper.get("grants", []),
            "referenced_works": paper.get("referenced_works", []),
            "ids": paper.get("ids", []),
            "counts_by_year": paper.get("counts_by_year", []),
        }
        for paper in response
    ]


def _works_generator(
    mailto: str,
    perpage: str,
    oa_id: Union[str, List[str]],
    filter_criteria: Union[str, List[str]],
    session: requests.Session,
) -> Iterator[list]:
    """Creates a generator that yields a list of works from the OpenAlex API based on a
    given work ID.

    Args:
        mailto (str): The email address to use for the API.
        perpage (str): The number of results to return per page.
        oa_id (Union[str, List[str]): The work ID to use for the API.
        filter_criteria (Union[str, List[str]]): The filter criteria to use for the API.
        session (requests.Session): The requests session to use.

    Yields:
        Iterator[list]: A generator that yields a list of works from the OpenAlex API
        based on a given work ID.
    """
    cursor = "*"
    assert isinstance(
        filter_criteria, type(oa_id)
    ), "filter_criteria and oa_id must be of the same type."

    # multiple filter criteria
    if isinstance(filter_criteria, list) and isinstance(oa_id, list):
        filter_string = ",".join(
            [f"{criteria}:{id_}" for criteria, id_ in zip(filter_criteria, oa_id)]
        )
    else:
        filter_string = f"{filter_criteria}:{oa_id}"

    cursor_url = (
        f"https://api.openalex.org/works?filter={filter_string}"
        f"&mailto={mailto}&per-page={perpage}&cursor={{}}"
    )

    try:
        # make a call to estimate total number of results
        response = session.get(cursor_url.format(cursor), timeout=20)
        data = response.json()

        while response.status_code == 429:  # needs testing (try with 200)
            logger.info("Waiting for 1 hour...")
            time.sleep(3600)
            response = session.get(cursor_url.format(cursor), timeout=20)
            data = response.json()

        logger.info("Fetching data for %s", oa_id[:50])
        total_results = data["meta"]["count"]
        num_calls = total_results // int(perpage) + 1
        logger.info("Total results: %s, in %s calls", total_results, num_calls)
        while cursor:
            response = session.get(cursor_url.format(cursor), timeout=20)
            data = response.json()
            results = data.get("results")
            cursor = data["meta"].get("next_cursor", False)
            yield results

    except Exception as e:  # pylint: disable=broad-except
        logger.error("Error fetching data for %s: %s", oa_id, e)
        yield []


def preprocess_oa_ids(
    oa_ids: Union[str, List[str], Dict[str, str]], group_oa_ids: bool
) -> List[str]:
    """Preprocesses oa_ids to ensure they are in the correct format."""
    if isinstance(oa_ids, str):
        oa_ids = [oa_ids]
    if isinstance(oa_ids, dict):
        oa_ids = list(oa_ids.values())
    if group_oa_ids:
        oa_ids = list(_chunk_oa_ids(oa_ids))
    return oa_ids


def _chunk_oa_ids(ids: List[str], chunk_size: int = 50) -> Generator[str, None, None]:
    """Yield successive chunk_size-sized chunks from ids."""
    for i in range(0, len(ids), chunk_size):
        yield "|".join(ids[i : i + chunk_size])


def fetch_papers_for_id(
    oa_id: Union[str, List[str]],
    mailto: str,
    perpage: str,
    filter_criteria: Union[str, List[str]],
) -> List[dict]:
    """Fetches all papers cited by a specific work ID."""
    assert isinstance(
        filter_criteria, type(oa_id)
    ), "filter_criteria and oa_id must be of the same type."
    papers_for_id = []
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=0.3)
    session.mount("https://", HTTPAdapter(max_retries=retries))
    for page, papers in enumerate(
        _works_generator(mailto, perpage, oa_id, filter_criteria, session)
    ):
        papers_for_id.extend(_parse_results(papers))
        logger.info(
            "Fetching page %s. Total papers collected: %s",
            page,
            len(papers_for_id),
        )

    return papers_for_id


def yield_papers_for_id(
    oa_id: Union[str, List[str]],
    mailto: str,
    perpage: str,
    filter_criteria: Union[str, List[str]],
) -> Generator[Dict[str, List[dict]], None, None]:
    """Fetches all papers cited by a specific work ID."""
    assert isinstance(
        filter_criteria, type(oa_id)
    ), "filter_criteria and oa_id must be of the same type."
    papers_collected = 0
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=0.3)
    session.mount("https://", HTTPAdapter(max_retries=retries))
    for page, papers in enumerate(
        _works_generator(mailto, perpage, oa_id, filter_criteria, session)
    ):
        results = _parse_results(papers)
        papers_collected += len(results)
        logger.info(
            "Fetching page %s. Total papers collected: %s",
            page,
            papers_collected,
        )

        yield {f"p{page}": results}


def fetch_papers_eager(
    processed_ids: Union[List[str], List[List[str]]],
    mailto: str,
    perpage: int,
    filter_criteria: Union[str, List[str]],
    slice_keys: bool,
) -> Dict[str, List[dict]]:
    """Fetches papers eagerly."""
    return {
        oa_id if not slice_keys else f"s{str(i)}": fetch_papers_for_id(
            oa_id=oa_id,
            mailto=mailto,
            perpage=perpage,
            filter_criteria=filter_criteria,
        )
        for i, oa_id in enumerate(processed_ids)
    }


def fetch_papers_lazy(
    processed_ids: Union[List[str], List[List[str]]],
    mailto: str,
    perpage: int,
    filter_criteria: Union[str, List[str]],
    slice_keys: bool,
) -> Dict[str, Callable]:
    """Fetches papers lazily."""
    return {
        (
            oa_id if not slice_keys else f"s{str(i)}"
        ): lambda oa_id=oa_id: fetch_papers_for_id(
            oa_id=oa_id,
            mailto=mailto,
            perpage=perpage,
            filter_criteria=filter_criteria,
        )
        for i, oa_id in enumerate(processed_ids)
    }


def fetch_papers_parallel(
    processed_ids: Union[List[str], List[List[str]]],
    mailto: str,
    perpage: int,
    filter_criteria: Union[str, List[str]],
    parallel_jobs: int = 4,
) -> Dict[str, List[Callable]]:
    """Fetches papers in parallel."""
    # slice oa_ids
    oa_id_chunks = [processed_ids[i : i + 80] for i in range(0, len(processed_ids), 80)]
    logger.info("Slicing data. Number of oa_id_chunks: %s", len(oa_id_chunks))
    return {
        f"s{str(i)}": lambda chunk=chunk: Parallel(n_jobs=parallel_jobs, verbose=10)(
            delayed(fetch_papers_for_id)(oa_id, mailto, perpage, filter_criteria)
            for oa_id in chunk
        )
        for i, chunk in enumerate(oa_id_chunks)
    }


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


def chunk_list(input_list: List[str], chunk_size: int) -> List[List[str]]:
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


def retrieve_oa_works_chunk(
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
