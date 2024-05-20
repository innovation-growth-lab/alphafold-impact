"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
from kedro.io import AbstractDataset
from ..data_analysis_descriptive_translational.nodes import (  # pylint: disable=relative-beyond-top-level
    get_entrez_ptype_pmid,
)
from joblib import Parallel, delayed
from Bio import Entrez


Entrez.email = "david.ampudia@nesta.org.uk"
logger = logging.getLogger(__name__)


def get_sb_lab_outputs(
    data_loaders: AbstractDataset,
    mapping_df: AbstractDataset,
):
    """
    Retrieves outputs from data loaders and performs data processing.

    Args:
        data_loaders (AbstractDataset): A dictionary-like object containing data loaders.
        mapping_df (AbstractDataset): A dataset containing mapping information.

    Returns:
        pd.DataFrame: The processed data.

    """
    outputs = []
    for i, loader in enumerate(data_loaders.values()):
        logger.info("Loading data batch %d / %d", i + 1, len(data_loaders))
        data_batch = loader()

        # drop "display_name", "authorships"
        data_batch = data_batch.drop(
            columns=["display_name", "authorships", "counts_by_year", "ids"]
        )

        # transform mesh, concepts, topics to be a list of the ids (ie first item in sublists)
        data_batch["mesh_terms"] = data_batch["mesh_terms"].apply(
            lambda x: [y[0] for y in x] if x is not None else []
        )
        data_batch["concepts"] = data_batch["concepts"].apply(
            lambda x: [y[0] for y in x] if x is not None else []
        )
        data_batch["topics"] = data_batch["topics"].apply(
            lambda x: [y[0] for y in x] if x is not None else []
        )

        outputs.append(data_batch)
    data = pd.concat(outputs)

    logger.info("Merging data with mapping information")
    mapping_df.drop_duplicates(subset="author", inplace=True)
    data = data.merge(
        mapping_df[["author", "seed"]], left_on="pi_id", right_on="author", how="left"
    )

    data["publication_date"] = pd.to_datetime(data["publication_date"])

    return data


def compute_publication_production(data: pd.DataFrame):
    """
    Compute monthly and yearly counts of publications for each 'pi_id' and 'seed'.

    Args:
        data (pd.DataFrame): The input DataFrame containing publication data.

    Returns:
        tuple: A tuple containing two DataFrames - monthly_publications and
            yearly_publications.
            - monthly_publications
            - yearly_publications
    """

    # Set 'publication_date' as the index
    data.set_index("publication_date", inplace=True)

    # Group by 'pi_id' and 'seed', then resample to get monthly and yearly counts of publications
    grouped = data.groupby(["pi_id", "seed"])
    monthly_publications = grouped.resample("M").size().reset_index(name="count")
    yearly_publications = grouped.resample("Y").size().reset_index(name="count")

    return monthly_publications, yearly_publications


def preprocess_for_event_study(
    data: pd.DataFrame, level0: pd.DataFrame
) -> pd.DataFrame:
    """
    Preprocesses the data for event study analysis.

    Args:
        data (pd.DataFrame): The input data to be preprocessed.
        level0 (pd.DataFrame): The level0 data used for merging.

    Returns:
        pd.DataFrame: The preprocessed data.

    """
    level0["publication_date"] = pd.to_datetime(level0["publication_date"])
    level0["parent_publication_date"] = pd.to_datetime(
        level0["parent_publication_date"]
    )

    logger.info("Merging data with level0 data")
    merged_data = pd.merge(
        data,
        level0[["id", "parent_id", "parent_publication_date"]],
        on="id",
        how="inner",
    )

    logger.info("Sorting data by 'parent_id' and 'parent_publication_date'")
    merged_data.sort_values(
        by=["parent_id", "parent_publication_date"],
        ascending=[True, True],
        inplace=True,
    )

    logger.info(
        "Computing time difference between 'publication_date' and 'parent_publication_date'"
    )
    parent_id_in_list = merged_data["parent_id"].isin(
        ["W3177828909", "W3211795435", "W3202105508"]
    )

    logger.info("Reordering data based on 'parent_id'")
    merged_data = pd.concat(
        [merged_data.loc[parent_id_in_list], merged_data.loc[~parent_id_in_list]]
    )

    logger.info("Dropping duplicates based on 'pi_id'")
    merged_data.drop_duplicates(subset="pi_id", keep="first", inplace=True)

    logger.info("Mapping 'parent_publication_date' to 'pi_id'")
    pi_id_date_dict = dict(
        zip(merged_data["pi_id"], merged_data["parent_publication_date"])
    )
    data["parent_publication_date"] = data["pi_id"].map(pi_id_date_dict)

    logger.info(
        "Computing time difference between 'publication_date' and 'parent_publication_date'"
    )
    data["time"] = (
        data["publication_date"] - data["parent_publication_date"]
    ) / pd.Timedelta(days=90)

    # Convert to integer
    final_data = data.copy()
    final_data.dropna(subset=["time"], inplace=True)
    final_data["time"] = final_data["time"].astype(int)

    return final_data


def get_event_study_outputs(data: pd.DataFrame) -> tuple:
    """
    Compute event study outputs based on the provided data and level0 data.

    Args:
        data (pd.DataFrame): The input data.
        level0 (pd.DataFrame): The level0 data.

    Returns:
        pd.DataFrame: Primary data.
        pd.DataFrame: The computed event study outputs.
        pd.DataFrame: The computed event study citations.


    """
    final_data_counts = (
        data.groupby(["pi_id", "time", "seed"]).size().reset_index(name="count")
    )
    final_data_citations = (
        data.groupby(["pi_id", "time", "seed"])["cited_by_count"]
        .sum()
        .reset_index(name="count")
    )

    return final_data_counts, final_data_citations


def get_event_study_strength(data, sc_data_af, sc_data_ct):
    """
    Compute event study strength based on the provided data and strong links data.

    Args:
        data (pd.DataFrame): The input data.
        sc_data_af (pd.DataFrame): The strong links data for AlphaFold.
        sc_data_ct (pd.DataFrame): The strong links data for counterfactual papers.

    Returns:
        pd.DataFrame: The computed event study outputs.
        pd.DataFrame: The computed event study citations.
    """
    sc_data = pd.concat([sc_data_af, sc_data_ct])

    # drop none in strength
    sc_data = sc_data.dropna(subset=["strength"])
    # drop duplicates
    sc_data = sc_data.drop_duplicates(subset=["id"], keep="first")

    data["strong"] = data["id"].isin(sc_data["id"])

    # for PI with at least a strong link, set all papers to strong
    pi_strong = data[data["strong"]]["pi_id"].unique()
    data.loc[data["pi_id"].isin(pi_strong), "strong"] = True

    final_data = data.copy()
    final_data.dropna(subset=["time"], inplace=True)
    final_data["time"] = final_data["time"].astype(int)

    final_data_counts = (
        final_data.groupby(["pi_id", "time", "seed", "strong"])
        .size()
        .reset_index(name="count")
    )
    # drop authors with seed other & strong True → these are authors not above the AF
    # threshold who nonetheless cited it, and above the threshold in other SB papers.
    final_data_counts = final_data_counts[
        ~((final_data_counts["seed"] == "other") & (final_data_counts["strong"]))
    ]

    final_data_citations = (
        final_data.groupby(["pi_id", "time", "seed", "strong"])["cited_by_count"]
        .sum()
        .reset_index(name="count")
    )
    final_data_citations = final_data_citations[
        ~((final_data_citations["seed"] == "other") & (final_data_citations["strong"]))
    ]

    return final_data_counts, final_data_citations


def get_event_study_pdb_submissions(
    data: pd.DataFrame, pdb_data: pd.DataFrame
) -> tuple:
    """
    Compute event study outputs based on the provided data and pdb_data.

    Args:
        data (pd.DataFrame): The input data.
        pdb_data (pd.DataFrame): The pdb data.

    Returns:
        pd.DataFrame: The computed event study outputs.
        pd.DataFrame: The computed event study citations.
    """

    logger.info("Merging data with pdb_data data")

    data = pd.merge(
        data,
        pdb_data[["id", "resolution", "R_free"]],
        on="id",
        how="inner",
    )

    final_data = data.copy()
    final_data.dropna(subset=["time"], inplace=True)
    final_data["time"] = final_data["time"].astype(int)

    final_data_counts = (
        final_data.groupby(["pi_id", "time", "seed", "resolution", "R_free"])
        .size()
        .reset_index(name="count")
    )
    final_data_citations = (
        final_data.groupby(["pi_id", "time", "seed", "resolution", "R_free"])[
            "cited_by_count"
        ]
        .sum()
        .reset_index(name="count")
    )
    return final_data_counts, final_data_citations


def get_event_study_predictive_outputs(data: pd.DataFrame) -> tuple:
    """
    Calculate the counts and citations of event study predictive outputs.

    Args:
        data (pd.DataFrame): The input DataFrame containing the data.

    Returns:
        tuple: A tuple containing two DataFrames:
            - final_data_counts: DataFrame with counts of event study predictive outputs.
            - final_data_citations: DataFrame with citations of event study predictive outputs.
    """

    data["protein_concept"] = data["concepts"].apply(
        lambda x: (
            any(element == "C18051474" or element == "C47701112" for element in x)
        )
    )

    final_data = data.copy()
    final_data.dropna(subset=["time"], inplace=True)
    final_data["time"] = final_data["time"].astype(int)

    final_data_counts = (
        final_data.groupby(["pi_id", "time", "seed", "protein_concept"])
        .size()
        .reset_index(name="count")
    )
    final_data_citations = (
        final_data.groupby(["pi_id", "time", "seed", "protein_concept"])[
            "cited_by_count"
        ]
        .sum()
        .reset_index(name="count")
    )
    return final_data_counts, final_data_citations


def get_event_study_cc(data: pd.DataFrame, icite_data: pd.DataFrame):
    """
    Get event study for clinical citations.

    Args:
        data (pd.DataFrame): The input data containing information about articles.
        icite_data (pd.DataFrame): The input data containing iCite information.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: A tuple containing three DataFrames:
            - data_pmid: The merged data with iCite information.
            - final_data_counts: The counts of clinical citations grouped by various attributes.
            - final_data_citations: The sum of cited_by_count grouped by various attributes.
    """

    # preprocess icite data
    icite_data["pmid"] = icite_data["pmid"].astype(str)

    logger.info("Merging chains with iCite data")

    # merge on 'pmid'
    data_pmid = data.merge(icite_data[["pmid", "cited_by_clin"]], how="left", on="pmid")
    data_pmid = data_pmid.drop_duplicates(subset=["id"])

    logger.info("Exploding cited_by_clin")
    data_pmid = data_pmid[data_pmid["cited_by_clin"].astype(str) != "nan"]

    # drop if cited_by_clin is None
    data_pmid = data_pmid[data_pmid["cited_by_clin"].notnull()]

    data_pmid["cited_by_clin"] = data_pmid["cited_by_clin"].apply(
        lambda x: x.split(" ")
    )
    data_pmid = data_pmid.explode("cited_by_clin")

    logger.info("Creating clinical article links")
    data_pmid["ca_link"] = data_pmid["cited_by_clin"].apply(
        lambda x: f"https://pubmed.ncbi.nlm.nih.gov/{x}"
    )

    # reindex
    data_pmid.reset_index(inplace=True, drop=True)

    # drop dup
    data_pmid.drop_duplicates(subset=["id", "cited_by_clin"], inplace=True)

    def apply_func(row):
        return get_entrez_ptype_pmid(row)
    
    results = Parallel(n_jobs=8, verbose=10)(
        delayed(apply_func)(row) for row in data_pmid["cited_by_clin"]
    )

    data_pmid[["ca_publication_type", "ca_publication_date"]] = pd.DataFrame(results)

    # data_pmid[["ca_publication_type", "ca_publication_date"]] = data_pmid[
    #     "cited_by_clin"
    # ].apply(get_entrez_ptype_pmid)

    logger.info("Collected clinical citation PMIDs")
    final_data = data_pmid.copy()
    final_data.dropna(subset=["time"], inplace=True)
    final_data["time"] = final_data["time"].astype(int)

    final_data_counts = (
        final_data.groupby(["pi_id", "time", "seed", "ca_publication_type"])
        .size()
        .reset_index(name="count")
    )
    final_data_citations = (
        final_data.groupby(["pi_id", "time", "seed", "ca_publication_type"])[
            "cited_by_count"
        ]
        .sum()
        .reset_index(name="count")
    )
    return data_pmid, final_data_counts, final_data_citations

def get_event_study_pc(
    data: pd.DataFrame, patent_data: pd.DataFrame
):
    """
    Get event study data for protein concepts.

    Args:
        data (pd.DataFrame): The input data containing information about events.
        patent_data (pd.DataFrame): The patent data to be merged with the input data.

    Returns:
        tuple: A tuple containing three elements:
            - data (pd.DataFrame): The merged data.
            - final_data_counts (pd.DataFrame): The final data grouped by various columns and the count of occurrences.
            - final_data_citations (pd.DataFrame): The final data grouped by various columns and the sum of citations.
    """
    
    # merge patent_data on data
    data = data.merge(
        patent_data, how="inner", left_on="doi", right_on="NPL Resolved External ID(s)"
    )

    final_data = data.copy()
    final_data.dropna(subset=["time"], inplace=True)
    final_data["time"] = final_data["time"].astype(int)

    final_data_counts = (
        final_data.groupby(["pi_id", "time", "seed", "protein_concept"])
        .size()
        .reset_index(name="count")
    )
    final_data_citations = ( # TODO This should be patent citations, not article 
        final_data.groupby(["pi_id", "time", "seed", "protein_concept"])[
            "cited_by_count"
        ]
        .sum()
        .reset_index(name="count")
    )

    return data, final_data_counts, final_data_citations