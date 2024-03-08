"""
This module contains functions for fetching and preprocessing data from
the OA API.

OpenAlex baseline:
    collect_papers: Collect papers based on the provided work IDs.
    load_oa_ids: Load the file corresponding to a particular oa_id in a PartitionedDataSet,
        extract all ids, and return these as a list.

OpenAlex - Gateway to Research:
    preprocess_publication_doi: Preprocess the Gateway to Research publication data to include
        doi values that are compatible with OA filter module.
    create_list_doi_inputs: Create a list of doi values from the Gateway to Research publication
        data.
    load_referenced_oa_ids: Load referenced work IDs from the dataset.

OpenAlex - Citation Depth:
    fetch_citation_depth: Iterates over an updating list of papers to process, yielding the
        response from collect papers for each paper in the list. As papers are collected, they are
        added to the set, while new, one-level deeper papers, are added to the list.
    create_network_graph: Creates the network graph from the edges, a list of tuples with the format
        (target, source).
    
 OpenAlex works for concepts:
    retrieve_oa_works_for_concepts_and_years: Collect papers based on
        the provided concept IDs and years.
"""

import logging
from typing import List, Dict, Sequence, Union, Callable, Tuple, Generator
from itertools import chain
import pandas as pd
from kedro.io import AbstractDataset
from .utils import (
    preprocess_oa_ids,
    fetch_papers_eager,
    fetch_papers_lazy,
    fetch_papers_parallel,
    retrieve_oa_works_chunk,
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
    concepts: bool = False,
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

    Returns:
        Union[Dict[str, Callable], Dict[str, List[dict]]]: A dictionary containing
        the collected papers.
            If eager_loading is True, the values are Callables that fetch the papers.
            If eager_loading is False, the values are Lists of dictionaries representing the
                fetched papers.
    """
    # preprocess oa_ids
    oa_ids = preprocess_oa_ids(oa_ids, group_oa_ids)

    # if concepts, simplify code and run direct fetch
    if concepts:
        return yield_papers_for_id(oa_ids, mailto, perpage, filter_criteria)

    # fetch papers for each oa_id
    if not parallelise:
        if eager_loading:
            return fetch_papers_eager(
                oa_ids, mailto, perpage, filter_criteria, slice_keys
            )
        return fetch_papers_lazy(oa_ids, mailto, perpage, filter_criteria, slice_keys)
    return fetch_papers_parallel(oa_ids, mailto, perpage, filter_criteria)


def load_oa_ids(oa_id: str, dataset: Sequence[AbstractDataset]) -> List[str]:
    """
    Loads the file corresponding to a particular oa_id in a PartitionedDataSet,
    extracts all ids, and returns these as a list.

    Args:
        oa_id (str): The oa_id to load the file for.
        dataset (PartitionedDataSet): The PartitionedDataSet containing the files.

    Returns:
        List[str]: A list of all ids extracted from the file.
    """
    data = dataset[oa_id]()
    ids = [paper["id"].replace("https://openalex.org/", "") for paper in data]
    logger.info("Loaded %s ids for %s", len(ids), oa_id)
    return ids


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
    concept_id_chunks = chunk_list(concept_ids, chunk_size)

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


def load_referenced_oa_ids(
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
    oa_ids = set()
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
                oa_id = work.get("id").replace("https://openalex.org/", "")
                doi = work.get("doi").replace("https://doi.org/", "")
                oa_doi_dict[oa_id] = {}
                oa_doi_dict[oa_id]["doi"] = doi
                oa_doi_dict[oa_id]["referenced_works"] = ref_works
                oa_ids.update(ref_works)

    oa_ids = list(oa_ids)
    logger.info("Work IDs loaded: %s", len(oa_ids))

    return oa_ids, oa_doi_dict


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
        concepts=True,
    )


def fetch_subfield_and_logic(
    oa_main_concept_ids: List[str],
    oa_and_concept_ids: List[List[str]],
    from_publication_date: str,
    api_config: Dict[str, str],
) -> Generator[Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]], None, None]:
    """
    Fetches papers based on a set of OA concepts that must at least have one concept
    in a separate list of concepts. Allows to control for AI research that is also
    related to fields of SB.

    Args:
        oa_main_concept_ids (List[str]): List of main concept IDs.
        oa_and_concept_ids (List[List[str]]): List of lists of concept IDs.
        from_publication_date (str): Publication date filter.
        api_config (Dict[str, str]): API configuration.

    Yields:
        Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]]: A generator that yields
            a tuple containing dictionaries of dataframes and dictionaries of lists
            of dictionaries.
    """

    # flatten the list of lists
    oa_main_concept_ids = "|".join(oa_main_concept_ids)
    oa_and_concept_ids = "|".join(list(chain(*oa_and_concept_ids)))

    # filter criteria
    filter_ids = [from_publication_date, oa_main_concept_ids, oa_and_concept_ids]

    return collect_papers(
        mailto=api_config["mailto"],
        perpage=api_config["perpage"],
        oa_ids=filter_ids,
        filter_criteria=["from_publication_date", "concepts.id", "concepts.id"],
        concepts=True,
    )
