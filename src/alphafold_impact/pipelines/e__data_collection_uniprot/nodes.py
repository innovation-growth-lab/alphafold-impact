"""
This is a boilerplate pipeline 'data_collection_pdb'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict, List
import requests
from requests.adapters import HTTPAdapter, Retry
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

logger = logging.getLogger(__name__)


def _fetch_pbd_uniprot_id(pdb_id) -> List[str]:
    """Fetch the list of PDB ids to be used for data collection"""
    # Fetch the list of PDB ids
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=0.1)
    session.mount("https://", HTTPAdapter(max_retries=retries))
    response = session.get(
        f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}", timeout=10
    )
    response = response.json()
    uniprot_obj = response[pdb_id.lower()].get("UniProt", {})

    # get the uniprot ids
    uniprot_ids = list(uniprot_obj.keys())

    # create a dataframe with the uniprot ids and the pdb id
    uniprot_ids = pd.DataFrame(uniprot_ids, columns=["uniprot_id"])
    uniprot_ids["pdb_id"] = pdb_id

    return uniprot_ids


def compute_complexity(data: Dict) -> float:
    """
    This function computes the complexity of a protein based on the following metrics:
    - Sequence Length
    - Sequence Variants
    - Secondary Structure
    - PTMs
    - GO Terms
    - PDB References

    The complexity is calculated as a weighted sum of the above metrics.

    Args:
        data (dict): JSON data from UniProt API

    Returns:
        float: Composite complexity score
    """
    # Example weights (adjust as needed based on importance)
    weights = {
        "sequence_length": 0.3,
        "sequence_variants": 0.1,
        "secondary_structure": 0.4,
        "ptms": 0.2,
        "go_terms": 0.2,
        "pdb_references": 0.3,
    }

    # Extract data from JSON
    complexity_metrics = {
        "sequence_length": data.get("sequence", {}).get("length", 0),
        "sequence_variants": data.get("extraAttributes", {})
        .get("countByFeatureType", {})
        .get("Sequence conflict", 0),
        "secondary_structure": (
            data.get("extraAttributes", {})
            .get("countByFeatureType", {})
            .get("Beta strand", 0)
            + data.get("extraAttributes", {})
            .get("countByFeatureType", {})
            .get("Helix", 0)
            + data.get("extraAttributes", {})
            .get("countByFeatureType", {})
            .get("Turn", 0)
        ),
        "ptms": (
            data.get("extraAttributes", {})
            .get("countByFeatureType", {})
            .get("Disulfide bond", 0)
            + data.get("extraAttributes", {})
            .get("countByFeatureType", {})
            .get("Glycosylation", 0)
            + data.get("extraAttributes", {})
            .get("countByFeatureType", {})
            .get("Lipidation", 0)
        ),
        "go_terms": len(
            [
                ref
                for ref in data.get("uniProtKBCrossReferences", [])
                if ref.get("database") == "GO"
            ]
        ),
        "pdb_references": sum(
            1
            for ref in data.get("uniProtKBCrossReferences", [])
            if ref.get("database") == "PDB"
        ),
    }

    # Composite Complexity Calculation
    composite_complexity = (
        weights["sequence_length"] * complexity_metrics["sequence_length"]
        + weights["sequence_variants"] * complexity_metrics["sequence_variants"]
        + weights["secondary_structure"] * complexity_metrics["secondary_structure"]
        + weights["ptms"] * complexity_metrics["ptms"]
        + weights["go_terms"] * complexity_metrics["go_terms"]
        + weights["pdb_references"] * complexity_metrics["pdb_references"]
    )

    return composite_complexity


def _extract_data(input_data):
    function_string = []
    triples_list = []

    for item in input_data:
        # collect all "FUNCTION" comment values
        if item.get("commentType") == "FUNCTION" and "texts" in item:
            for text_item in item["texts"]:
                if "value" in text_item:
                    function_string.append(text_item["value"])

        # extract triples (id, acronym, MIM ID) for "DISEASE" type comments
        if item.get("commentType") == "DISEASE" and "disease" in item:
            disease = item["disease"]
            disease_id = disease.get("diseaseId")
            acronym = disease.get("acronym")
            mim_id = disease.get("diseaseCrossReference", {}).get("id")
            if disease_id and acronym and mim_id:
                triples_list.append((disease_id, acronym, mim_id))

    # join the collected FUNCTION values into a single string
    joined_function_string = " ".join(function_string)

    return joined_function_string, triples_list


def safe_get(d, keys, default=np.nan):
    """
    Safely retrieves a nested key from a dictionary, returning a
        default value if any key is missing.
    """
    for key in keys:
        if isinstance(d, dict):
            d = d.get(key, None)
        else:
            return default
    return d if d is not None else default


def get_uniprot_details(uniprot_id: str) -> Dict:
    """
    Get the details of a uniprot id from the uniprot API and return a dictionary
    that contains the following information:

    - uniprot_id
    - organism_id
    - score
    - protein_existence
    - protein_name
    - complexity

    Args:
        uniprot_id (str): Uniprot ID

    Returns:
        Dict: Dictionary containing the details of the uniprot id
    """
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=0.1)
    session.mount("https://", HTTPAdapter(max_retries=retries))
    response = session.get(
        f"https://rest.uniprot.org/uniprotkb/search?query=accession_id:{uniprot_id}",
        timeout=10,
    )
    response = response.json()["results"][0]

    function_string, disease_triples = _extract_data(response.get("comments", []))

    obj_ret = {
        "uniprot_id": uniprot_id,
        "organism_id": response.get("organism", {}).get("taxonId", np.nan),
        "score": response.get("annotationScore", np.nan),
        "protein_existence": response.get("proteinExistence", np.nan),
        "protein_name": (
            response.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", np.nan)
        ),
        "complexity": compute_complexity(response),
        "function_string": function_string,
        "disease_triples": disease_triples,
    }

    return obj_ret


def fetch_pdb_uniprot_map(pdb_data: pd.DataFrame) -> pd.DataFrame:
    """
    This function fetches the mapping between PDB and UniProt IDs.

    Args:
        pdb_data (pd.DataFrame): The DataFrame containing PDB details.

    Returns:
        pd.DataFrame: The DataFrame containing the mapping between PDB
          and UniProt IDs.
    """
    uniprot_ids = Parallel(n_jobs=6, verbose=10)(
        delayed(_fetch_pbd_uniprot_id)(pdb_id) for pdb_id in pdb_data["rcsb_id"]
    )
    uniprot_ids = pd.concat(uniprot_ids)

    return uniprot_ids


def fetch_uniprot_data(uniprot_ids: pd.DataFrame) -> pd.DataFrame:
    """
    Fetch the details of the uniprot ids from the uniprot API.

    Args:
        uniprot_ids (pd.DataFrame): The DataFrame containing the uniprot ids.

    Returns:
        pd.DataFrame: The DataFrame containing the details of the uniprot ids.
    """
    # get the unique ids to iterate requests for
    unique_uniprot_ids = uniprot_ids["uniprot_id"].unique()

    logger.info("Fetching details for %s uniprot ids", len(unique_uniprot_ids))

    # paralleling fetching uniprot details
    uniprot_details = Parallel(n_jobs=6, verbose=10)(
        delayed(get_uniprot_details)(uniprot_id) for uniprot_id in unique_uniprot_ids
    )

    # create a dataframe
    uniprot_details = pd.DataFrame(uniprot_details)

    # merge the data
    uniprot_data = uniprot_ids.merge(uniprot_details, on="uniprot_id", how="left")

    return uniprot_data
