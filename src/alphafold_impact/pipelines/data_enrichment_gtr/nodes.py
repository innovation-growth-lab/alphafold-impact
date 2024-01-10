"""
This module contains the functions used to enrich the Gateway to Research
(GtR) data with OpenAlex data.

Functions:
    create_pipeline: Creates the data enrichment pipeline.
    preprocess_publication_doi: Preprocess the Gateway to Research publication data to include
        doi values that are compatible with OA filter module.
    create_list_doi_inputs: Create a list of doi values from the Gateway to Research publication data.
    create_dictionary_doi_to_oa: Create a dictionary of doi values to OpenAlex publication data.
"""
from typing import Dict, List, Any
import pandas as pd


def preprocess_publication_doi(df: pd.DataFrame) -> pd.DataFrame:
    """Preprocess the Gateway to Research publication data to include
    doi values that are compatible with OA filter module.

    Args:
        df (pd.DataFrame): The Gateway to Research publication data.

    Returns:
        pd.DataFrame: The preprocessed publication data.
    """
    df["doi"] = (
        df["doi"]
        .str.replace("/dx.", "/", regex=True)
        .str.replace("http:", "https:", regex=False)
    )
    return df


def create_list_doi_inputs(df: pd.DataFrame) -> list:
    """Create a list of doi values from the Gateway to Research publication data.

    Args:
        df (pd.DataFrame): The Gateway to Research publication data.

    Returns:
        list: A list of doi values.
    """
    return df[df["doi"].notnull()]["doi"].tolist()


def create_dictionary_doi_to_oa(
    dictionary: Dict[str, List[Dict[str, Any]]]
) -> Dict[str, str]:
    """Create a dictionary of doi values to OpenAlex publication data.

    Args:
        dictionary (Dict[str, List[Dict[str, Any]]]): The OpenAlex publication data.

    Returns:
        Dict[str, str]: A dictionary of doi values to OpenAlex publication data.
    """
    doi_to_oa_dict = {}
    for _, value in dictionary.items():
        for item in value:
            doi_to_oa_dict[item["doi"]] = item["id"].replace(
                "https://openalex.org/", ""
            )
    return doi_to_oa_dict
