"""
This module contains functions for fetching and preprocessing data from the OA API.

OpenAlex baseline:
    collect_papers: Collect papers based on the provided work IDs.
    load_work_ids: Load the file corresponding to a particular work_id in a PartitionedDataSet,
        extract all ids, and return these as a list.

OpenAlex - Gateway to Research:
    preprocess_publication_doi: Preprocess the Gateway to Research publication data to include
        doi values that are compatible with OA filter module.
    create_list_doi_inputs: Create a list of doi values from the Gateway to Research publication
        data.
    load_referenced_work_ids: Load referenced work IDs from the dataset.

OpenAlex - Citation Depth:
    fetch_citation_depth: Iterates over an updating list of papers to process, yielding the
        response from collect papers for each paper in the list. As papers are collected, they are
        added to the set, while new, one-level deeper papers, are added to the list.
    create_network_graph: Creates the network graph from the edges, a list of tuples with the format
        (target, source).
"""

import logging
from typing import List, Dict, Sequence, Union, Callable, Tuple
import pandas as pd
import networkx as nx
from kedro.io import AbstractDataset
from joblib import Parallel, delayed
from .utils import (
    preprocess_work_ids,
    fetch_papers_eager,
    fetch_papers_lazy,
    fetch_papers_parallel,
    retrieve_oa_works_chunk,
)

logger = logging.getLogger(__name__)


def collect_papers(
    mailto: str,
    perpage: str,
    work_ids: Union[str, List[str], Dict[str, str]],
    filter_criteria: str,
    group_work_ids: bool = False,
    eager_loading: bool = False,
    slice_keys: bool = False,
    parallelise: bool = False,
) -> Union[Dict[str, Callable], Dict[str, List[dict]]]:
    """
    Collects papers based on the provided work IDs.

    Args:
        mailto (str): The email address to be used for API requests.
        perpage (str): The number of papers to fetch per page.
        work_ids (Union[str, List[str], Dict[str, str]]): The work IDs to collect papers for.
        filter_criteria (str): The filter to apply when fetching papers.
        group_work_ids (bool, optional): Whether to group the work IDs. Defaults to False.
        eager_loading (bool, optional): Whether to eagerly load all papers. Defaults to False.
        slice_keys (bool, optional): Whether to use slices as keys in the result dictionary.
            Defaults to False.
        parallelise (bool, optional): Whether to parallelize the fetching of papers.
            Defaults to False.

    Returns:
        Union[Dict[str, Callable], Dict[str, List[dict]]]: A dictionary containing
        the collected papers.
            If eager_loading is True, the values are Callables that fetch the papers.
            If eager_loading is False, the values are Lists of dictionaries representing the
                fetched papers.
    """
    # preprocess work_ids
    work_ids = preprocess_work_ids(work_ids, group_work_ids)

    # fetch papers for each work_id
    if not parallelise:
        if eager_loading:
            return fetch_papers_eager(
                work_ids, mailto, perpage, filter_criteria, slice_keys
            )
        return fetch_papers_lazy(work_ids, mailto, perpage, filter_criteria, slice_keys)
    return fetch_papers_parallel(work_ids, mailto, perpage, filter_criteria)


def load_work_ids(work_id: str, dataset: Sequence[AbstractDataset]) -> List[str]:
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
        chunk_size (int, optional): The number of concept IDs to process in each API call
            (default is 40).
        per_page (int, optional): The number of results to retrieve per API call (default is 200).

    Returns:
        pd.DataFrame: A pandas DataFrame containing all unique retrieved works.
    """
    all_works = []
    concept_id_chunks = _chunk_list(concept_ids, chunk_size)

    for index, concept_chunk in enumerate(concept_id_chunks, start=1):
        logger.info("Processing chunk %s/%s", index, len(concept_id_chunks))
        chunk_works = retrieve_oa_works_chunk(
            concept_chunk, publication_years, per_page
        )
        all_works.extend(chunk_works)

    all_works_df = pd.DataFrame(all_works).drop_duplicates(subset=["id"])
    logger.info("Retrieved %s works.", len(all_works_df))
    return all_works_df


def preprocess_publication_doi(df: pd.DataFrame) -> pd.DataFrame:
    """Preprocess the Gateway to Research publication data to include
    doi values that are compatible with OA filter module.

    Args:
        df (pd.DataFrame): The Gateway to Research publication data.

    Returns:
        pd.DataFrame: The preprocessed publication data.
    """
    if "doi" in df.columns:
        df["doi"] = (
            df["doi"]
            .str.replace("/dx.", "/", regex=True)
            .str.replace("http:", "https:", regex=False)
            .str.split()
            .str[0]
            .apply(lambda x: x if isinstance(x, str) and x.count(":") <= 1 else None)
        )
    return df


def create_list_doi_inputs(df: pd.DataFrame) -> list:
    """Create a list of doi values from the Gateway to Research publication data.

    Args:
        df (pd.DataFrame): The Gateway to Research publication data.

    Returns:
        list: A list of doi values.
    """
    return df[df["doi"].notnull()]["doi"].drop_duplicates().tolist()


def load_referenced_work_ids(
    dataset: Dict[str, Sequence[Sequence[AbstractDataset]]]
) -> Tuple[List[str], Dict[str, str]]:
    """
    Load referenced work IDs from the dataset.

    Args:
        dataset (Dict[str, Sequence[Sequence[AbstractDataset]]]): A dictionary
            containing the dataset.

    Returns:
        Tuple[List[str], Dict[str, str]]: A tuple containing the list of work IDs
            and a dictionary mapping work IDs to DOIs.
    """
    oa_doi_dict = {}
    work_ids = set()
    for timestamp, loader in dataset.items():
        logger.info("Loading work IDs from %s", timestamp)
        data = loader()

        # for each OR call that we made to openalex
        for call in data:
            for work in call:
                ref_works = [
                    w.replace("https://openalex.org/", "")
                    for w in work["referenced_works"]
                ]
                work_id = work.get("id").replace("https://openalex.org/", "")
                doi = work.get("doi").replace("https://doi.org/", "")
                oa_doi_dict[work_id] = {}
                oa_doi_dict[work_id]["doi"] = doi
                oa_doi_dict[work_id]["referenced_works"] = ref_works
                work_ids.update(ref_works)

    work_ids = list(work_ids)
    logger.info("Work IDs loaded: %s", len(work_ids))

    return work_ids, oa_doi_dict


def fetch_citation_depth(
    seed_paper: str, api_config: Dict[str, str], filter_config: str
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]]:
    """
    Iterates over an updating list of papers to process, yielding the response from collect
    papers for each paper in the list. As papers are collected, they are added to the set, while
    new, one-level deeper papers, are added to the list.

    Args:
        seed_paper (str): The seed paper to start from.
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
        child_papers = Parallel(n_jobs=8, backend="loky", verbose=10)(
            delayed(collect_papers)(
                work_ids=paper,
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


def create_network_graph(edges: pd.DataFrame) -> nx.Graph:
    """
    Creates the network graph from the edges, a list of tuples with the format (target, source).

    Args:
        edges (pd.DataFrame): The edges of the network graph.

    Returns:
        nx.Graph: The network graph.
    """
    G = nx.from_pandas_edgelist(edges, "target", "source")
    return G
