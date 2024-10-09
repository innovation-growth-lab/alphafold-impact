"""
This module contains functions for fetching and preprocessing data from
the OpenAlex (OA) API.

Functions:
    collect_papers(mailto, perpage, oa_ids, filter_criteria, group_oa_ids=False, **kwargs):
        Collects papers based on the provided work IDs.
    fetch_subfield_baseline(oa_concept_ids, from_publication_date, api_config):
        Fetches the baseline data for a subfield based on the given concept IDs and 
        publication date.

This module also contains functions for fetching citation data to a specific depth and enriching
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

Functions:
    fetch_citation_to_specific_depth(seed_paper, api_config, filter_config, max_depth, ...):
        Fetches citations to a specific depth from a seed paper.
    preprocess_baseline_data(data, alphafold_papers):
        Preprocesses the baseline data by extracting the paper IDs from the given dataset and 
        removing the IDs that are present in the `alphafold_papers` list.
    fetch_level_2_ct_papers(ct_data):
        Fetches the level 2 parent CT papers from a DataFrame.
    _process_flatten_dict(child_papers_dict):
        Processes and flattens a dictionary of child papers, removing entries before 2021 and 
        cleaning IDs.
"""

import logging
from typing import List, Dict, Union, Callable, Tuple, Generator
import pandas as pd
from joblib import Parallel, delayed
from kedro.io import AbstractDataset
from .utils import (
    preprocess_oa_ids,
    fetch_papers_eager,
    fetch_papers_lazy,
    fetch_papers_parallel,
    yield_papers_for_id,
)

logger = logging.getLogger(__name__)


def collect_papers(
    mailto: str,
    perpage: str,
    oa_ids: Union[str, List[str], List[List[str]], Dict[str, str]],
    filter_criteria: Union[str, List[str]],
    group_oa_ids: bool = False,
    eager_loading: bool = False,
    slice_keys: bool = False,
    parallelise: bool = False,
    bypass: bool = False,
    skip_preprocess: bool = False,
    **kwargs,
) -> Union[Dict[str, Callable], Dict[str, List[dict]]]:
    """
    Collects papers based on the provided work IDs.

    Args:
        mailto (str): The email address to be used for API requests.
        perpage (str): The number of papers to fetch per page.
        oa_ids (Union[str, List[str], List[List[str]], Dict[str, str]]): The work IDs to fetch
            papers for.
        filter_criteria (Union[str, List[str]]): The filter to apply when fetching papers.
        group_oa_ids (bool, optional): Whether to group the work IDs. Defaults to False.
        eager_loading (bool, optional): Whether to eagerly load all papers. Defaults to False.
        slice_keys (bool, optional): Whether to use slices as keys in the result dictionary.
            Defaults to False.
        parallelise (bool, optional): Whether to parallelize the fetching of papers.
            Defaults to False.
        bypass (bool, optional): Whether the work IDs are concepts. Defaults to False.

    Returns:
        Union[Dict[str, Callable], Dict[str, List[dict]]]: A dictionary containing
        the collected papers.
            If eager_loading is True, the values are Callables that fetch the papers.
            If eager_loading is False, the values are Lists of dictionaries representing the
                fetched papers.
    """
    if not skip_preprocess:
        # preprocess oa_ids
        oa_ids = preprocess_oa_ids(oa_ids, group_oa_ids)

    # if concepts, simplify code and run direct fetch
    if bypass:
        return yield_papers_for_id(oa_ids, mailto, perpage, filter_criteria)

    # fetch papers for each oa_id
    if not parallelise:
        if eager_loading:
            return fetch_papers_eager(
                oa_ids, mailto, perpage, filter_criteria, slice_keys, **kwargs
            )
        return fetch_papers_lazy(oa_ids, mailto, perpage, filter_criteria, slice_keys)
    return fetch_papers_parallel(oa_ids, mailto, perpage, filter_criteria)


def fetch_subfield_baseline(
    oa_concept_ids: List[str],
    from_publication_date: str,
    api_config: Dict[str, str],
) -> Generator[Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]], None, None]:
    """
    Fetches the baseline data for a subfield based on the given concept IDs and publication date.

    Args:
        oa_concept_ids (List[str]): List of concept IDs.
        from_publication_date (str): The publication date from which to fetch the papers.
        api_config (Dict[str, str]): API configuration parameters.

    Yields:
        Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]]: A generator that yields a tuple
            containing two dictionaries.
            The first dictionary contains dataframes with different fields as keys.
            The second dictionary contains lists of dictionaries with concept information.
    """

    # Preprocess concept IDs
    oa_concept_ids = "|".join(oa_concept_ids)

    # filter criteria
    filter_ids = [from_publication_date, oa_concept_ids]

    # fetch papers
    return collect_papers(
        mailto=api_config["mailto"],
        perpage=api_config["perpage"],
        oa_ids=filter_ids,
        filter_criteria=["from_publication_date", "concepts.id"],
        bypass=True,
    )


def fetch_citation_to_specific_depth(
    seed_paper: Union[str, List[str]],
    api_config: Dict[str, str],
    filter_config: str,
    max_depth: int,
    start_level: int = 0,
    papers_seen: set = None,
):
    """
    Fetches citations to a specific depth from a seed paper.

    Args:
        seed_paper (Union[str, List[str]]): The seed paper to start fetching
            citations from.
        api_config (Dict[str, str]): API configuration parameters.
        filter_config (str): Filter configuration for fetching papers.
        max_depth (int): The maximum depth to fetch citations to.

    Yields:
        dict: A dictionary containing the fetched papers at each level of depth.
    """
    if papers_seen is None:
        papers_seen = set()
    processed_papers, papers_to_process, level = (
        papers_seen,
        set(seed_paper),
        start_level,
    )

    while level < max_depth:
        next_level_papers = set()
        batch_num = 0

        while papers_to_process:
            batch_num += 1
            current_batch = {
                papers_to_process.pop() for _ in range(min(800, len(papers_to_process)))
            }
            logger.info("Processing %d papers", len(current_batch))

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
            yield {f"{level}/p{batch_num}": new_papers}

        papers_to_process = next_level_papers
        level += 1


def preprocess_baseline_data(
    data: AbstractDataset, alphafold_papers: List[str]
) -> List[str]:
    """
    Preprocesses the baseline data by extracting the paper IDs from the given dataset
    and removing the IDs that are present in the `alphafold_papers` list.

    Args:
        data (AbstractDataset): The dataset containing the JSON data.
        alphafold_papers (List[str]): The list of paper IDs to be excluded.

    Returns:
        List[str]: The list of paper IDs after preprocessing.

    """
    seed_papers = set()
    for loader in data.values():
        json_data = loader()
        for record in json_data:
            paper = record.get("id", "").replace("https://openalex.org/", "")
            seed_papers.add(paper)
    seed_papers -= set(alphafold_papers)
    return list(seed_papers)


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
