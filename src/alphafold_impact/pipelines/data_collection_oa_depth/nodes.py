"""This module contains functions for fetching citation data to a specific depth and enriching
the data with additional MeSH tags.

The enrichment process is done in the following steps:
    1. Load a batch of data from the input dataset.
    2. Load JSON data, transform it into a DataFrame, and perform data transformations.
    3. Fetch additional mesh terms for DOIs in a DataFrame using Biopython's Entrez module.
        - Filter DataFrame for rows where mesh_terms is empty.
        - Create list of DOIs.
        - Fetch PubMed IDs for each DOI.
        - Use PubMed IDs to fetch mesh terms from full-data objects.
        - Update the DataFrame with the new mesh terms for the corresponding DOI, 
            removing the [DOI] part.

OpenAlex - Citation Depth:
    fetch_citation_all_depth: Iterates over an updating list of papers to process, yielding the
        response from collect papers for each paper in the list. As papers are collected, they are
        added to the set, while new, one-level deeper papers, are added to the list.
    fetch_citation_specific_depth: Iterates over an updating list of papers to process, yielding
        the response from collect papers for each paper in the list. As papers are collected, they
        are added to the set, while new, one-level deeper papers, are added to the list. Now
        includes level tracking to limit depth of citation tree exploration.
    _process_flatten_dict: Process and flatten, removing entries before 2021 and cleaning IDs.

OpenAlex - Biopython MeSH tagging:
    data_generator: Load a batch of data from the input dataset.
    json_loader: Load JSON data, transform it into a DataFrame, and perform data transformations.
    fetch_additional_mesh: Fetches additional mesh terms for DOIs in a DataFrame.
"""

import logging
from typing import Dict, List, Tuple, Generator
import pandas as pd
from joblib import Parallel, delayed
from ..data_collection_oa.nodes import collect_papers  # pylint: disable=E0402

logger = logging.getLogger(__name__)


def fetch_citation_all_depth(
    seed_paper: str, api_config: Dict[str, str], filter_config: str
) -> Generator[Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]], None, None]:
    """
    Iterates over an updating list of papers to process, yielding the response from collect
    papers for each paper in the list. As papers are collected, they are added to the set, while
    new, one-level deeper papers, are added to the list.

    Args:
        seed_paper (str): The seed work ID paper to start from, ie. AlphaFold's.
        api_config (Dict[str, str]): The API configuration.
        filter_config (str): The filter to apply when fetching papers.

    Yields:
        Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]]: A tuple containing the edges and
        the papers.
    """
    processed_paper_ids = set()
    papers_to_process = {seed_paper}  # Use a set for uniqueness and efficient look-up

    while papers_to_process:
        logger.info("Processing %s papers", len(papers_to_process))
        current_batch = set()

        # take up to 200 papers for processing
        while papers_to_process and len(current_batch) < 200:
            current_batch.add(papers_to_process.pop())

        processed_paper_ids.update(current_batch)
        logger.info("Processing the following papers: %s", current_batch)

        # parallel fetching of papers
        child_papers = Parallel(n_jobs=6, backend="loky", verbose=10)(
            delayed(collect_papers)(
                oa_ids=paper,
                mailto=api_config["mailto"],
                perpage=api_config["perpage"],
                filter_criteria=filter_config,
                eager_loading=True,
            )
            for paper in current_batch
        )

        # flatten the list of dicts into a single dict
        child_papers_flat = {k: v for d in child_papers for k, v in d.items()}

        lengths = {key: len(value) for key, value in child_papers_flat.items()}
        logger.info("Lengths of value lists: %s", lengths)

        new_papers = set()
        for parent, children in child_papers_flat.items():
            edge_list = []

            # removing papers published before 2021-01-01
            children = [
                child
                for child in children
                if child.get("publication_date") >= "2021-01-01"
            ]

            # adding edges
            edge_list = [
                (parent, child.get("id", "").replace("https://openalex.org/", ""))
                for child in children
            ]

            # updating the list of new papers
            new_papers.update(
                [
                    clean_id
                    for child in children
                    if (
                        clean_id := child.get("id", "").replace(
                            "https://openalex.org/", ""
                        )
                    )
                    not in processed_paper_ids
                ]
            )

            edge_list_df = pd.DataFrame(edge_list, columns=["target", "source"])
            yield {parent: edge_list_df}, {parent: children}

        # extend the papers_to_process without duplicates
        papers_to_process.update(new_papers - processed_paper_ids)


def fetch_citation_to_specific_depth(
    seed_paper: str, api_config: Dict[str, str], filter_config: str, max_depth: int
):
    """
    Fetches citations to a specific depth from a seed paper.

    Args:
        seed_paper (str): The seed paper to start fetching citations from.
        api_config (Dict[str, str]): API configuration parameters.
        filter_config (str): Filter configuration for fetching papers.
        max_depth (int): The maximum depth to fetch citations to.

    Yields:
        dict: A dictionary containing the fetched papers at each level of depth.
    """
    processed_papers, papers_to_process, level = set(), {seed_paper}, 0

    while level < max_depth:
        next_level_papers = set()
        # outputs = {}

        while papers_to_process:
            current_batch = {
                papers_to_process.pop() for _ in range(min(800, len(papers_to_process)))
            }
            logger.info("Processing %d papers", len(current_batch))

            child_papers = Parallel(n_jobs=8, backend="loky", verbose=10)(
                delayed(collect_papers)(
                    oa_ids=paper,
                    mailto=api_config["mailto"],
                    perpage=api_config["perpage"],
                    filter_criteria=filter_config,
                    eager_loading=True,
                )
                for paper in current_batch
            )

            processed_papers.update(current_batch)
            child_papers_flat = {k: v for d in child_papers for k, v in d.items()}
            new_papers = _process_flatten_dict(child_papers_flat)

            next_level_papers.update(
                source["id"]
                for sources in new_papers.values()
                for source in sources
                if source["id"] not in processed_papers
            )
            # outputs.update(new_papers)
            yield {f"{level}/{key}": value for key, value in new_papers.items()}

        papers_to_process = next_level_papers
        level += 1


def _process_flatten_dict(child_papers_dict: Dict[str, any]) -> Dict[str, any]:
    """Process and flatten, removing entries before 2021 and cleaning IDs."""
    return {
        target: [
            {
                **{k: v for k, v in source.items() if k != "referenced_works"},
                "id": source.get("id", "").replace("https://openalex.org/", ""),
            }
            for source in sources
            if source.get("publication_date", "0") >= "2021-01-01"
        ]
        for target, sources in child_papers_dict.items()
        if any(
            source.get("publication_date", "0") >= "2021-01-01" for source in sources
        )
    }
