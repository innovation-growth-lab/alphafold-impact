"""
This module contains the functions used to enrich the Gateway to Research
(GtR) data with OpenAlex data.

Functions:
    preprocess_publication_doi: Preprocess the Gateway to Research publication data to include
        doi values that are compatible with OA filter module.
    create_list_doi_inputs: Create a list of doi values from the Gateway to Research publication data.
    test_node: Test node.
    process_for_partition_storing: Process the OpenAlex publication data for storing in partitions.
"""
import logging
from typing import Dict, List, Sequence, Tuple
import pandas as pd
from kedro.io import AbstractDataset

logger = logging.getLogger(__name__)


def preprocess_publication_doi(df: pd.DataFrame) -> pd.DataFrame:
    """Preprocess the Gateway to Research publication data to include
    doi values that are compatible with OA filter module.

    Args:
        df (pd.DataFrame): The Gateway to Research publication data.

    Returns:
        pd.DataFrame: The preprocessed publication data.
    """
    # Ensure the DOI column exists
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
    return df[(df["doi"].notnull()) & (df.duplicated(subset="doi"))]["doi"].tolist()[:50]


def load_referenced_work_ids(
    dataset: Dict[str, Sequence[Sequence[AbstractDataset]]]
) -> Tuple[List[str], Dict[str, str]]:
    """
    Load referenced work IDs from the dataset.

    Args:
        dataset (Dict[str, Sequence[Sequence[AbstractDataset]]]): A dictionary containing the dataset.

    Returns:
        Tuple[List[str], Dict[str, str]]: A tuple containing the list of work IDs and a dictionary mapping work IDs to DOIs.
    """
    oa_doi_dict = {}
    work_ids = set()
    for timestamp, loader in dataset.items():
        logger.info("Loading work IDs from %s", timestamp)
        data = loader()

        # for each OR call that we made to openalex
        for call in data:
            for work in call:
                work_id = work.get("id").replace("https://openalex.org/", "")
                doi = work.get("doi").replace("https://doi.org/", "")
                oa_doi_dict[work_id] = doi
                ref_works = [
                    w.replace("https://openalex.org/", "")
                    for w in work["referenced_works"]
                ]
                work_ids.update(ref_works)

    work_ids = list(work_ids)
    logger.info("Work IDs loaded: %s", len(work_ids))

    return work_ids, oa_doi_dict


def test_node(data_inputs):
    logger.info("Testing")
    return data_inputs[list(data_inputs.keys())[-1]]
