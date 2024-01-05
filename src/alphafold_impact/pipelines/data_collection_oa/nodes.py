"""
This module contains functions for fetching and preprocessing data from the GtR API.

Functions:
- _revert_abstract_index(abstract_inverted_index: Dict[str, Sequence[int]]) -> str: 
    Reverts the abstract inverted index to the original text.
- _parse_results(response: List[Dict]) -> Dict[str, List[str]]: 
    Parses OpenAlex API response to retain specific information.
- _citation_works_generator(mailto: str, perpage: str, work_id: str, direction: str) -> Iterator[list]: 
    Creates a generator that yields a list of works from the OpenAlex API based on a given work ID.
- collect_citation_papers(mailto: str, perpage: str, work_id: str, direction: str) -> list: 
    Collects all papers cited by specific work IDs.
"""

import logging
from typing import Iterator, List, Dict, Sequence, Union
from kedro.io import AbstractVersionedDataset
from requests.adapters import HTTPAdapter, Retry
import requests
from typing import List

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
        work_ids (List[str]): A list of work IDs to filter by either 'cites' or 'cited_by'.
        direction (str): The direction of the citation. Either 'cites' or 'cited_by'.

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


def load_work_ids(work_id: str, dataset: Sequence[AbstractVersionedDataset]) -> List[str]:
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