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
import random
from typing import Iterator, List, Dict, Sequence, Union, Callable, Generator
import time
from requests.adapters import HTTPAdapter, Retry
import requests
from joblib import Parallel, delayed

logger = logging.getLogger(__name__)


MAIL_TO_CANDIDATES = [
    "david.ampudia@nesta.org.uk",
    "data_analytics@nesta.org.uk",
    "david.ampudia@bse.eu",
    "david.ampudia@upf.edu",
    "yanyan.leung@nesta.org.uk",
    "layla.gemmell@nesta.org.uk",
    "george.richardson@nesta.org.uk",
    "hugo.cuello@nesta.org.uk",
    "nyangala.zolho@nesta.org.uk",
    "edoardo.trimarchi@nesta.org.uk"
]


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
            "id": paper.get("id", "").replace("https://openalex.org/", ""),
            "doi": paper.get("doi", ""),
            "publication_date": paper.get("publication_date", ""),
            "authorships": paper.get("authorships", []),
            "cited_by_count": paper.get("cited_by_count", ""),
            "concepts": paper.get("concepts", []),
            "mesh_terms": paper.get("mesh", []),
            "topics": paper.get("topics", []),
            "grants": paper.get("grants", []),
            "ids": paper.get("ids", []),
            "counts_by_year": paper.get("counts_by_year", []),
            "citation_normalized_percentile": paper.get(
                "citation_normalized_percentile", []
            ),
            "cited_by_percentile_year": paper.get("cited_by_percentile_year", []),
            "fwci": paper.get("fwci", ""),
        }
        for paper in response
    ]


def _make_request(session: requests.Session, url: str) -> dict:
    response = session.get(url, timeout=20)
    while response.status_code == 429:
        logger.info("Waiting for 1 hour...")
        time.sleep(30)  # 1 hour wait
        response = session.get(url, timeout=20)
    return response.json()


def _log_and_yield_results(
    session: requests.Session, url: str, num_calls: int
) -> Iterator[list]:
    for page in range(1, num_calls + 1):
        data = _make_request(session, url.format(page))
        results = data.get("results")
        yield results


def _works_generator(
    perpage: str,
    oa_id: Union[str, List[str]],
    filter_criteria: Union[str, List[str]],
    session: requests.Session,
    sample_size: int = -1,
) -> Iterator[list]:
    """Creates a generator that yields a list of works from the OpenAlex API based on a
    given work ID.

    Args:
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

    if isinstance(filter_criteria, list) and isinstance(oa_id, list):
        filter_string = ",".join(
            [f"{criteria}:{id_}" for criteria, id_ in zip(filter_criteria, oa_id)]
        )
    else:
        filter_string = f"{filter_criteria}:{oa_id}"

    mailto = random.choice(MAIL_TO_CANDIDATES)

    if sample_size == -1:
        cursor_url = (
            f"https://api.openalex.org/works?filter={filter_string}"
            f"&mailto={mailto}&per-page={perpage}&cursor={{}}"
        )
        try:
            data = _make_request(session, cursor_url.format(cursor))
            total_results = data["meta"]["count"]
            num_calls = total_results // int(perpage) + 1
            logger.info("Total results: %s, in %s calls", total_results, num_calls)
            while cursor:
                data = _make_request(session, cursor_url.format(cursor))
                results = data.get("results")
                cursor = data["meta"].get("next_cursor", False)
                yield results
        except Exception as e: # pylint: disable=broad-except
            logger.error("Error fetching data for %s: %s", oa_id, e)
            yield []
    else:
        cursor_url = (
            f"https://api.openalex.org/works?filter={filter_string}&seed=123"
            f"&mailto={mailto}&per-page={perpage}&sample={sample_size}&page={{}}"
        )
        try:
            data = _make_request(session, cursor_url.format(1))
            total_results = data["meta"]["count"]
            num_calls = total_results // int(perpage) + 1
            logger.info("Total results: %s, in %s calls", total_results, num_calls)
            yield from _log_and_yield_results(session, cursor_url, num_calls)
        except Exception as e: # pylint: disable=broad-except
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
    perpage: str,
    filter_criteria: Union[str, List[str]],
    **kwargs,
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
        _works_generator(perpage, oa_id, filter_criteria, session, **kwargs)
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
        _works_generator(perpage, oa_id, filter_criteria, session)
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
    perpage: int,
    filter_criteria: Union[str, List[str]],
    slice_keys: bool,
    **kwargs,
) -> Dict[str, List[dict]]:
    """Fetches papers eagerly."""
    if len(processed_ids) == len(filter_criteria):
        processed_ids = [processed_ids]
    return {
        oa_id if not slice_keys else f"s{str(i)}": fetch_papers_for_id(
            oa_id=oa_id,
            perpage=perpage,
            filter_criteria=filter_criteria,
            **kwargs,
        )
        for i, oa_id in enumerate(processed_ids)
    }


def fetch_papers_lazy(
    processed_ids: Union[List[str], List[List[str]]],
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
            perpage=perpage,
            filter_criteria=filter_criteria,
        )
        for i, oa_id in enumerate(processed_ids)
    }


def fetch_papers_parallel(
    processed_ids: Union[List[str], List[List[str]]],
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
            delayed(fetch_papers_for_id)(oa_id, perpage, filter_criteria)
            for oa_id in chunk
        )
        for i, chunk in enumerate(oa_id_chunks)
    }
