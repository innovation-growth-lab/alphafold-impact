"""Pipeline nodes for the data_collection_s2 pipeline.

Functions:
    load_alphafold_citation_ids: Loads the citation IDs for the AlphaFold paper.
    fetch_citation_details: Fetches citation details for a given work ID.
    get_citation_details: Retrieves citation details for a given list of work IDs 
        and filters them based on a specified AlphaFold DOI.
"""

import logging
from typing import Sequence, Dict, Tuple, Union
import re
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from kedro.io import AbstractDataset

logger = logging.getLogger(__name__)


def load_alphafold_citation_ids(
    input_loaders: AbstractDataset, work_id: str
) -> Sequence[Tuple[str]]:
    """Loads the citation IDs for the AlphaFold paper. Requires OA pipeline."""

    logger.info("Loading citation IDs for Alphafold")
    citations = input_loaders[work_id]()

    # iterate over four possible ids (oa, doi, mag, pmid)
    ids = []
    for citation in citations:
        id_dict = citation.get("ids", {})
        openalex = id_dict.get("openalex", "").replace("https://openalex.org/", "")
        doi = re.sub(r".*?(10\..*)", "\\1", id_dict.get("doi", ""))
        mag = id_dict.get("mag", "")
        pmid = (
            id_dict.get("pmid", "").split("/")[-1]
            if "/" in id_dict.get("pmid", "")
            else id_dict.get("pmid", "")
        )
        ids.append((openalex, doi, mag, pmid))

    return ids


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


def get_alphafold_citation_details(
    work_ids: Sequence[Tuple[str, str, str, str]], af_doi: str, **kwargs
) -> pd.DataFrame:
    """Retrieves citation details for a given list of work IDs and filters
    them based on the specified AlphaFold DOI.

    Args:
        work_ids (Sequence[Tuple[str, str, str, str]]): A list of work IDs.
        af_doi (str): The AlphaFold DOI to filter the citation details.
        **kwargs: Additional keyword arguments to be passed to the
            fetch_citation_details function.

    Returns:
        pd.DataFrame: A DataFrame containing the citation details for the
        specified work IDs, filtered by the AlphaFold DOI. The DataFrame
        has the following columns:
        - oa: The work ID.
        - citedPaper: The details of the cited paper, including the title,
          authors, and publication information.
        - otherColumns: Any additional columns present in the citation details.

    """
    citation_details = {}
    for work_id_tuple in work_ids:
        oa, doi, mag, pmid = work_id_tuple
        logger.info("Fetching citation details for %s", oa)
        for prefix, id_ in [("DOI:", doi), ("PMID:", pmid), ("MAG:", mag)]:
            if not id_:
                logger.warning("No relevant %s id", prefix[:-1])
                continue
            work_id = f"{prefix}{id_}"
            try:
                data = fetch_citation_details(work_id, **kwargs)
                for item in data:
                    if (eids := item["citedPaper"]["externalIds"]) is not None:
                        if eids.get("DOI", "") == af_doi:
                            logger.info("Found citation details using %s", work_id)
                            citation_details[oa] = {
                                "item": item,
                                "doi": doi,
                                "mag": mag,
                                "pmid": pmid,
                            }
                            break
                else:
                    logger.warning("No relevant citation found in %s", work_id)
                    continue
                break
            except Exception:  # pylint: disable=broad-except
                logger.warning("Failed to fetch citation details for %s", work_id)

    # explode multiple citations
    citation_details = process_citations(citations=citation_details, parent=af_doi)
    return {af_doi: citation_details}, [tuple([af_doi, 34265844])]


#  list(set(citation_details["child_doi"].to_list()))


def fetch_citation_strength(
    parent_data: Sequence[Dict[str, AbstractDataset]],
    logs: Sequence[Tuple[str, str]],
    **kwargs,
) -> pd.DataFrame:
    """
    Fetches the citation strength for a new level based on the parent data and logs.

    Args:
        parent_data (Sequence[Dict[str, AbstractDataset]]): The parent data containing 
            the datasets.
        logs (Sequence[Tuple[str, str]]): The logs containing the paper IDs.
        **kwargs: Additional keyword arguments.

    Returns:
        pd.DataFrame: The citation strength data.

    Yields:
        Dict[str, List[Dict[str, Any]]], List[List[str]]: The citation details for each 
            paper and the processed paper IDs.
    """
    logger.info("Fetching citation strength for a new level")
    processed_paper_ids = set(tuple(log) for log in logs)
    paper_ids_to_process = set()

    for _, loader in parent_data.items():
        data = loader()
        set_of_papers = (
            set(zip(data["child_doi"].to_list(), data["child_pmid"].to_list()))
            - processed_paper_ids
        )
        paper_ids_to_process.update(set_of_papers)

    while paper_ids_to_process:
        current_paper = paper_ids_to_process.pop()
        logger.info(
            "Fetching citation details for %s. Papers remaining: %s",
            current_paper,
            len(paper_ids_to_process),
        )
        try:
            citation_details = fetch_citation_details(
                work_id=f"DOI:{current_paper[0]}", **kwargs
            )
        except requests.exceptions.HTTPError:
            pass

        if not citation_details and current_paper[1]:
            try:
                citation_details = fetch_citation_details(
                    work_id=f"PMID:{current_paper[1]}", **kwargs
                )
            except requests.exceptions.HTTPError:
                continue

        logger.info("Fetched citation details for %s", current_paper)
        citation_details = process_citations(
            citations=citation_details, parent=current_paper
        )

        processed_paper_ids.add(current_paper)

        logger.info("Fetched citation details for %s", current_paper)
        yield {current_paper[0]: citation_details}, [current_paper]

    logger.info(
        "Finished fetching citation strength for a new level. Flushing IDs list."
    )

    # back to list of lists
    processed_paper_ids = [list(log) for log in processed_paper_ids]

    yield pd.DataFrame(), processed_paper_ids


def process_citations(
    citations: Union[Dict[str, Union[str, Sequence[str]]], Sequence],
    parent: Union[str, Tuple[str, str]],
) -> pd.DataFrame:
    """
    Process the citations data and convert it into a pandas DataFrame.

    Args:
        citations (Sequence[Dict[str, Union[str, Sequence[str]]]]): A sequence of
            dictionaries representing the citations data.
            Each dictionary should have the following keys:
            - 'intents': A list of strings representing the intents.
            - 'contexts': A list of strings representing the contexts.
        parent (str): Parent paper that receives the citation.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the processed citations data.
            The DataFrame has the following columns:
            - 'parent': The parent identifier.
            - 'child_oa': The child identifier.
            - 'child_doi': The DOI of the child paper.
            - 'child_mag': The MAG ID of the child paper.
            - 'child_pmid': The PMID of the child paper.+
            - 'influential': Whether the citation is influential or not.
            - 'intent': The intent of the citation.
            - 'context': The context of the citation.
    """
    rows = []
    if isinstance(citations, dict):
        for child_oa, details in citations.items():
            nested_details = details.get("item", {})
            context_intents = nested_details.get("contextsWithIntent", [])
            influential = nested_details.get("isInfluential", "")
            child_doi = details.get("doi", "")
            child_mag = details.get("mag", "")
            child_pmid = details.get("pmid", "")
            for item in context_intents:
                context = item.get("context", "")
                intents = item.get("intents") or [""]
                for intent in intents:
                    rows.append(
                        [
                            parent,
                            "",
                            child_oa,
                            child_doi,
                            child_mag,
                            child_pmid,
                            "",
                            influential,
                            intent,
                            context,
                        ]
                    )

    elif isinstance(citations, list):
        for citation in citations:
            parent_doi = parent[0]
            parent_pmid = parent[1]
            context_intents = citation.get("contextsWithIntent", [])
            influential = citation.get("isInfluential", "")
            ids = citation.get("citingPaper", {}).get("externalIds", {})
            if not ids:
                # skip if no ids
                continue
            child_doi = ids.get("DOI", "")
            child_mag = ids.get("MAG", "")
            child_pmid = ids.get("PubMed", "")
            child_cid = ids.get("CorpusId", "")
            for item in context_intents:
                context = item.get("context", "")
                intents = item.get("intents") or [""]
                for intent in intents:
                    rows.append(
                        [
                            parent_doi,
                            parent_pmid,
                            "",
                            child_doi,
                            child_mag,
                            child_pmid,
                            child_cid,
                            influential,
                            intent,
                            context,
                        ]
                    )

    # return dataframe
    return pd.DataFrame(
        rows,
        columns=[
            "parent_doi",
            "parent_pmid",
            "child_oa",
            "child_doi",
            "child_mag",
            "child_pmid",
            "child_cid",
            "influential",
            "intent",
            "context",
        ],
    )
