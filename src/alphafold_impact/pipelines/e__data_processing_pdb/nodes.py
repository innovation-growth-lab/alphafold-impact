"""
This module contains the functions to process PDB data.

The main functions in this module are:

* collect_pdb_details - Collects PDB details from a DataFrame and API configuration.
* merge_uniprot_data - Merges Uniprot data with PDB data.
* process_similarity_data - Process the entire similarity data file in chunks and 
    compute novelty metrics.
"""

import logging
import pandas as pd
import numpy as np
from .utils import (
    get_papers,
    preprocess_pdb_dates,
    filter_and_compute_metrics,
)

logger = logging.getLogger(__name__)


def collect_pdb_details(pdb_df: pd.DataFrame, api_config: dict) -> pd.DataFrame:
    """
    Collects PDB details from a DataFrame and API configuration.

    Args:
        pdb_df (pd.DataFrame): The DataFrame containing PDB details.
        api_config (dict): The API configuration.

    Returns:
        pd.DataFrame: The combined DataFrame with collected PDB details.
    """

    # get the set of unique pmid
    pmids = (
        pdb_df["pmid"]
        .replace(["", "None", "nan"], np.nan)
        .dropna()
        .astype(float)
        .astype(int)
        .astype(str)
        .unique()
        .tolist()
    )

    processed_papers_with_pmids = get_papers(pmids, api_config, label="ids.pmid")

    logger.info(
        "Processed %d/%d papers with pmids from openalex",
        len(processed_papers_with_pmids),
        len(pmids),
    )

    # filter out from pdb_df the rows with pmid in processed_papers_with_pmids
    spdb_df = pdb_df[
        ~pdb_df["pmid"].isin(
            processed_papers_with_pmids["pmid"].astype(float).astype(str)
        )
    ]

    # turns out, missing dois still have a pdb entry doi,
    # not included in the API response!
    spdb_df["doi"] = spdb_df.apply(
        lambda row: (
            f"10.2210/pdb{row['rcsb_id'].lower()}/pdb"
            if row["doi"] in ["", "None", "nan"]
            else row["doi"]
        ),
        axis=1,
    )

    logger.info(
        "Attempting to process %d papers with dois from openalex",
        len(spdb_df),
    )

    # extract the list of doi from spdb_df
    dois = spdb_df["doi"].unique().tolist()

    processed_papers_with_dois = get_papers(dois, api_config, label="doi")

    logger.info(
        "Processed %d/%d papers with dois from openalex",
        len(processed_papers_with_dois),
        len(dois),
    )

    # concatenate the two frames, drop duplicate ids
    oa_pdb_df = pd.concat(
        [processed_papers_with_pmids, processed_papers_with_dois]
    ).drop_duplicates(subset=["id"])

    pdb_with_pmid = pdb_df[~pdb_df["pmid"].isin(["", "None", "nan"])]
    pdb_with_pmid["pmid"] = pdb_with_pmid["pmid"].astype(float).astype(int).astype(str)

    logger.info("Merging openalex data with pdb data")
    oa_pdb_df_pmid = oa_pdb_df.merge(
        pdb_with_pmid[["rcsb_id", "pmid", "resolution", "R_free"]],
        how="right",
        on="pmid",
    )

    # remove rows with missing ids
    oa_pdb_df_pmid.dropna(subset="id", inplace=True)

    def _is_pdb_doi(doi):
        return doi is not None and doi.startswith("10.2210/pdb")

    # remove duplicate rows if one is the pdb entry
    oa_pdb_df_pmid["is_pdb_doi"] = oa_pdb_df_pmid["doi"].apply(_is_pdb_doi)
    oa_pdb_df_pmid = (
        oa_pdb_df_pmid.sort_values(by="is_pdb_doi")
        .drop_duplicates(subset="rcsb_id", keep="first")
        .drop(columns=["is_pdb_doi"])
    )

    logger.info("Merging openalex data with pdb data")
    oa_pdb_df_doi = oa_pdb_df.merge(
        pdb_df[["rcsb_id", "doi", "resolution", "R_free"]], how="right", on="doi"
    )

    logger.info("Concatenating openalex data with pdb data")
    combined_df = pd.concat([oa_pdb_df_pmid, oa_pdb_df_doi])  # .drop_duplicates("id")
    combined_df.dropna(subset="id", inplace=True)

    # remove again duplicate rows if one is the pdb entry
    combined_df["is_pdb_doi"] = combined_df["doi"].apply(_is_pdb_doi)
    combined_df = (
        combined_df.sort_values(by="is_pdb_doi")
        .drop_duplicates(subset="rcsb_id", keep="first")
        .drop(columns=["is_pdb_doi"])
    )

    logger.info("Finished processing %d papers", len(combined_df))
    # replace openalex.org in id
    combined_df["id"] = combined_df["id"].str.replace("https://openalex.org/", "")

    # drop ids
    combined_df.drop(columns=["ids"], inplace=True)

    return combined_df



def aggregate_to_pdb_level(similarity_chunks: pd.DataFrame) -> pd.DataFrame:
    """
    Distill similarity data to PDB-to-PDB level by extracting PDB IDs from
        query and target.

    Args:
        similarity_df (pd.DataFrame): The DataFrame containing similarity data.

    Returns:
        pd.DataFrame: The DataFrame containing computed metrics.
    """
    results = []
    for i, chunk in enumerate(similarity_chunks):
        logger.info("Processing chunk %d", i + 1)
        chunk["query"] = chunk["query"].str[:4]
        chunk["target"] = chunk["target"].str[:4]

        # Group by PDB-to-PDB comparisons and compute metrics
        pdb_level_chunk = chunk.groupby(["query", "target"], as_index=False).agg(
            {"alntmscore": "mean"}
        )

        results.append(pdb_level_chunk)

    # Concatenate all individual results
    pdb_level_df = pd.concat(results, ignore_index=True)

    # [HACK] Group again in case chunks split a query in 2+
    pdb_level_df = pdb_level_df.groupby(
        ["query", "target"], as_index=False
    ).agg(
        {
            "alntmscore": "mean"
        }
    )


    return pdb_level_df


def process_similarity_data(
    similarity_df: pd.DataFrame,
    pdb_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Process the entire similarity data file in chunks and compute novelty metrics.

    Args:
        similarity_file_path (str): The path to the similarity data file.
        pdb_df (pd.DataFrame): The DataFrame containing PDB details.
        chunksize (int): The size of each chunk to load.

    Returns:
        pd.DataFrame: The DataFrame containing computed metrics.
    """
    pdb_df = pdb_df[["id", "rcsb_id", "publication_date", "doi", "R_free", "resolution"]]

    logger.info("Preprocessing PDB dates")
    pdb_dates = preprocess_pdb_dates(pdb_df)

    logger.info("Filtering and computing metrics")
    results = filter_and_compute_metrics(similarity_df, pdb_dates)

    logger.info("Merging results with PDB data")
    pdb_df["rcsb_id"] = pdb_df["rcsb_id"].str.upper()
    pdb_df = pdb_df.merge(results, how="left", left_on="rcsb_id", right_on="query")

    return pdb_df


def merge_uniprot_data(pdb_df: pd.DataFrame, uniprot_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merges Uniprot data with PDB data.

    Args:
        pdb_df (pd.DataFrame): The DataFrame containing PDB details.
        uniprot_df (pd.DataFrame): The DataFrame containing Uniprot details.

    Returns:
        pd.DataFrame: The combined DataFrame with PDB and Uniprot details.
    """

    logger.info("Transformations to PDB data")
    pdb_df["R_free"] = pd.to_numeric(pdb_df["R_free"], errors="coerce")
    pdb_df["resolution"] = pd.to_numeric(pdb_df["resolution"], errors="coerce")

    pdb_df = pdb_df.assign(
        publication_date=pd.to_datetime(pdb_df["publication_date"]),
        quarter=pd.PeriodIndex(pdb_df["publication_date"], freq="Q"),
    )

    pdb_oa_means = pdb_df.groupby("id").agg(
        resolution_mean=("resolution", "mean"),
        R_free_mean=("R_free", "mean"),
        mean_tmscore=("mean_tmscore", "mean"),
        max_tmscore=("max_tmscore", "mean"),
        normalised_mean_tmscore=("normalised_mean_tmscore", "mean"),
        normalised_max_tmscore=("normalised_max_tmscore", "mean"),
    )

    logger.info("merging to uniprot to get per-PDB primary status")
    merged_uniprot_df = uniprot_df.merge(
        pdb_df, how="left", left_on="pdb_id", right_on="rcsb_id"
    )

    # create column primary_submission for a given uniprot_id
    merged_uniprot_df["primary_submission"] = (
        merged_uniprot_df.groupby("uniprot_id")["publication_date"]
        .transform("min")
        .eq(merged_uniprot_df["publication_date"])
    )

    # create unique assignment
    primary_pdbs = (
        merged_uniprot_df[["rcsb_id", "primary_submission"]]
        .sort_values("primary_submission")
        .drop_duplicates(subset=["rcsb_id"], keep="last")
    )

    # transform boolean to int
    primary_pdbs["primary_submission"] = primary_pdbs["primary_submission"].astype(
        "Int64"
    )

    # merge back
    pdb_df = pdb_df.merge(primary_pdbs, how="left")

    # merge uniprot data with pdb data
    intermediate_df = uniprot_df.merge(
        pdb_df, how="left", left_on="pdb_id", right_on="rcsb_id"
    )

    logger.info("Creating disease counts in intermediate_df")
    intermediate_df["num_diseases"] = intermediate_df["disease_triples"].apply(
        lambda x: len(x) if x is not None else 0
    )

    logger.info("Creating inverse frequency of organism occurrence")
    organism_frequencies = (
        merged_uniprot_df.groupby(["organism_id", "quarter"])["organism_id"]
        .nunique()
        .groupby(level=0)
        .cumsum()
        .reset_index(name="cumulative_unique_organism_count")
    )
    organism_frequencies["organism_rarity"] = (
        1 / organism_frequencies["cumulative_unique_organism_count"]
    )

    # merge organism_rarity to intermediate_df
    intermediate_df = intermediate_df.merge(
        organism_frequencies[["organism_id", "quarter", "organism_rarity"]],
        how="left",
        on=["organism_id", "quarter"],
    )

    logger.info("Creating publication-level UniProt features")
    oa_structural_df = (
        intermediate_df.groupby("id")
        .agg(
            num_uniprot_structures=("uniprot_id", "nunique"),
            num_pdb_ids=("rcsb_id", "nunique"),
            num_primary_submissions=("primary_submission", "sum"),
            score_mean=("score", "mean"),
            complexity_sum=("complexity", "sum"),
            complexity_mean=("complexity", "mean"),
            organism_rarity_mean=("organism_rarity", "mean"),
            organism_rarity_max=("organism_rarity", "max"),
            num_diseases=("num_diseases", "sum"),
        )
        .reset_index()
    )

    # merge grouped_df with pdb_oa_means
    oa_structural_df = oa_structural_df.merge(
        pdb_oa_means, how="left", left_on="id", right_index=True
    )

    return oa_structural_df
