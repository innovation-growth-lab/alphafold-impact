"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
from Bio import Entrez
from ..data_processing_chains.nodes import (  # pylint: disable=E0402
    get_entrez_ptype_pmid,
)

Entrez.email = "david.ampudia@nesta.org.uk"
logger = logging.getLogger(__name__)


def load_input_data(
    data: pd.DataFrame,
    source: str,
):
    """
    Load input data and perform necessary transformations.

    Args:
        data (pd.DataFrame): The input data to be loaded.
        source (str): The source of the data.

    Returns:
        pd.DataFrame: The loaded data with additional columns.

    """
    data["level"] = data["level"].astype(str)

    logger.info("Filter for level 0")
    data = data[data["level"] == "0"]

    logger.info("Create columns with source")
    data["source"] = source

    return data


def merge_inputs(**kwargs):
    """Merge inputs."""
    merged_data = pd.concat(
        [kwargs["alphafold_data"], kwargs["ct_data"], kwargs["other_data"]]
    )

    merged_data["dataset"] = merged_data.groupby("id")["source"].transform(",".join)

    # remove repeated "af", "ct", or "other" statements in dataset
    merged_data["dataset"] = merged_data["dataset"].apply(
        lambda x: ",".join(sorted(list(set(x.split(",")))))
    )

    merged_data = merged_data.drop_duplicates(subset="id", keep="first")

    return merged_data


def preprocess_sb_data(data: pd.DataFrame):
    """Preprocess the data."""
    # create a dictionary mapping id to publication_date
    id_date_dict = data.set_index("id")["publication_date"].to_dict()

    # create a new column 'parent_publication_date' by mapping the 'parent_id' column to the dictionary
    data["parent_publication_date"] = data["parent_id"].map(id_date_dict)

    # if parent_id is W3177828909, set parent_publication_date to 2021-06-25
    data.loc[data["parent_id"] == "W3177828909", "parent_publication_date"] = (
        "2021-06-25"
    )

    # change level to -1 for parent_id "W3177828909", "W3211795435", "W3202105508"
    data.loc[
        data["id"].isin(["W3177828909", "W3211795435", "W3202105508"]), "level"
    ] = -1

    return data


def get_cc_counts(data: pd.DataFrame, icite_data: pd.DataFrame):
    """
    Get citation counts for the given data.

    Args:
        data (pd.DataFrame): The input data containing information about articles.
        icite_data (pd.DataFrame): The input data containing citation information.

    Returns:
        pd.DataFrame: The combined data with citation counts.

    """
    # preprocess icite data
    icite_data["pmid"] = icite_data["pmid"].astype(str)

    logger.info("Merging chains with iCite data")

    # merge on 'pmid'
    data_pmid = data.merge(icite_data[["pmid", "cited_by_clin"]], how="left", on="pmid")
    data_pmid = data_pmid.drop_duplicates(
        subset=["parent_id", "id", "level", "dataset"]
    )

    # merge on 'doi'
    data_doi = data.merge(icite_data[["doi", "cited_by_clin"]], how="left", on="doi")
    data_doi = data_doi.drop_duplicates(subset=["parent_id", "id", "level", "dataset"])

    combined_data = pd.concat([data_pmid, data_doi]).drop_duplicates(
        subset=["parent_id", "id", "level", "dataset"]
    )

    logger.info("Exploding cited_by_clin")
    combined_data = combined_data[combined_data["cited_by_clin"].astype(str) != "nan"]

    # drop if cited_by_clin is None
    combined_data = combined_data[combined_data["cited_by_clin"].notnull()]

    combined_data["cited_by_clin"] = combined_data["cited_by_clin"].apply(
        lambda x: x.split(" ")
    )
    combined_data = combined_data.explode("cited_by_clin")

    logger.info("Creating clinical article links")
    combined_data["ca_link"] = combined_data["cited_by_clin"].apply(
        lambda x: f"https://pubmed.ncbi.nlm.nih.gov/{x}"
    )

    # reindex
    combined_data.reset_index(inplace=True, drop=True)

    # drop dup
    combined_data.drop_duplicates(subset=["id", "cited_by_clin"], inplace=True)

    combined_data["publication_type"] = combined_data["cited_by_clin"].apply(
        get_entrez_ptype_pmid
    )

    return combined_data


# strong and weak
# how many strong articles
# patenst
