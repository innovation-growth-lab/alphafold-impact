"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
from Bio import Entrez


Entrez.email = "david.ampudia@nesta.org.uk"
logger = logging.getLogger(__name__)


def get_entrez_ptype_pmid(pmid: str) -> str:
    """
    Retrieves the publication type for a given PubMed ID (PMID).

    Args:
        pmid (str): The PubMed ID (PMID) of the publication.

    Returns:
        str: The publication type of the specified publication.
    """
    try:
        stream = Entrez.efetch(db="pubmed", id=pmid, retmax="1")
        record = Entrez.read(stream)
        pub_type = str(
            record["PubmedArticle"][0]["MedlineCitation"]["Article"].get(
                "PublicationTypeList"
            )[0]
        )

        year = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][
            0
        ].get("Year")
        month = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][
            0
        ].get("Month")
        day = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][
            0
        ].get("Day")

        pub_date = f"{year}-{month}-{day}"
        return pd.Series([pub_type, pub_date])

    except Exception:  # pylint: disable=broad-except
        return pd.Series([None, None])


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
    data = data[(data["level"] == "0") | (data["level"] == "-1")]

    logger.info("Create columns with source")
    data["source"] = source

    return data

def load_input_applied_data(
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

    logger.info("Filter for level 1 and 2")
    data = data[(data["level"] == "0") |(data["level"] == "1") | (data["level"] == "2")]

    data = data[["id", "doi", "pmid", "parent_id", "publication_date", "parent_level", "level"]]

    # create a dictionary mapping id to publication_date
    id_date_dict = data.set_index("id")["publication_date"].to_dict()

    # create a new column 'parent_publication_date' by mapping the 'parent_id' column to the dictionary
    data["parent_publication_date"] = data["parent_id"].map(id_date_dict)

    # force level and parent_level to be str
    data["level"] = data["level"].astype(str)
    data["parent_level"] = data["parent_level"].astype(str)

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

    # create a row that is af's
    new_row = {
        "id": "W3177828909",
        "doi": "10.1038/s41586-021-03819-2",
        "pmid": "34265844",
        "title": "Highly accurate protein structure prediction with AlphaFold",
        "publication_date": "2021-07-15",
        "level": "-1",
        "parent_level": "-2",
        "source": "af",
    }

    data.loc[len(data)] = new_row

    # change level to -1 for parent_id "W3177828909", "W3211795435", "W3202105508"
    data.loc[data["id"].isin(["W3211795435", "W3202105508"]), "level"] = "-1"
    data.loc[data["id"].isin(["W3211795435", "W3202105508"]), "parent_level"] = "-2"

    # create a dictionary mapping id to publication_date
    id_date_dict = data.set_index("id")["publication_date"].to_dict()

    # create a new column 'parent_publication_date' by mapping the 'parent_id' column to the dictionary
    data["parent_publication_date"] = data["parent_id"].map(id_date_dict)

    # force level and parent_level to be str
    data["level"] = data["level"].astype(str)
    data["parent_level"] = data["parent_level"].astype(str)
    return data


def get_cc_papers(data: pd.DataFrame, icite_data: pd.DataFrame):
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

    combined_data[["ca_publication_type", "ca_publication_date"]] = combined_data[
        "cited_by_clin"
    ].apply(get_entrez_ptype_pmid)

    return combined_data


def get_cc_moderators(
    cc_data: pd.DataFrame,
    sc_data_af: pd.DataFrame,
    sc_data_ct: pd.DataFrame,
    pdb_data: pd.DataFrame,
):
    # Combine the 'id' columns of the three DataFrames into one Series
    combined_ids = pd.concat([sc_data_af["id"], sc_data_ct["id"]])

    # Check if the 'id' values in cc_data are in combined_ids
    cc_data["id_in_sc_data"] = cc_data["id"].isin(combined_ids)

    # Check if the 'id' values in cc_data are in pdb_data
    cc_data["id_in_pdb_data"] = cc_data["id"].isin(pdb_data["id"])

    # add R_free and resolution metrics
    cc_data = pd.merge(
        cc_data,
        pdb_data[["id", "R_free", "resolution"]],
        on="id",
        how="left",
    )


    # check if cited_by_clin is in pdb_data
    cc_data["cited_by_clin_in_pdb_data"] = cc_data["cited_by_clin"].isin(
        pdb_data["pmid"]
    )

    return cc_data


def get_patent_papers(data: pd.DataFrame, patent_data: pd.DataFrame):
    """
    Get patent counts for the given data.

    Args:
        data (pd.DataFrame): The input data containing information about articles.
        patent_data (pd.DataFrame): The input data containing patent information.

    Returns:
        pd.DataFrame: The combined data with patent counts.

    """
    # merge patent_data on data
    matched_data = data.merge(
        patent_data, how="inner", left_on="doi", right_on="NPL Resolved External ID(s)"
    )
    return matched_data


def get_patent_moderators(
    pc_data: pd.DataFrame,
    sc_data_af: pd.DataFrame,
    sc_data_ct: pd.DataFrame,
    pdb_data: pd.DataFrame,
):
    # Combine the 'id' columns of the three DataFrames into one Series
    combined_ids = pd.concat([sc_data_af["id"], sc_data_ct["id"]])

    # Check if the 'id' values in pc_data are in combined_ids
    pc_data["id_in_sc_data"] = pc_data["id"].isin(combined_ids)

    # add R_free and resolution metrics
    pc_data = pd.merge(
        pc_data,
        pdb_data[["id", "R_free", "resolution"]],
        on="id",
        how="left",
    )

    # Check if the 'id' values in pc_data are in pdb_data
    pc_data["id_in_pdb_data"] = pc_data["id"].isin(pdb_data["id"])

    return pc_data
