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
from ..data_analysis_lab_productivity.nodes import (  # pylint: disable=E0402
    # _preprocess_for_staggered_design,
    _get_pdb_activity,
    _get_patent_citations,
    _get_awards,
    _get_cc,
    _get_quarterly_aggregate_outputs,
    _get_time,
    _get_yearly_citations,
    _calculate_mesh_balance,
    _calculate_topic_share,
    _collect_covid_references
)

from Bio import Entrez

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

    logger.info("Filter for level 1 and 2")
    data = data[
        (data["level"] == "0") | (data["level"] == "1") | (data["level"] == "2")
    ]

    data = data[["id", "parent_id", "publication_date", "parent_level", "level"]]

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


def get_applied_lab_outputs(
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
        if i > 2:
            break
    data = pd.concat(outputs)

    logger.info("Merging data with mapping information")
    mapping_df.drop_duplicates(subset="author", inplace=True)
    data = data.merge(
        mapping_df[["author", "seed"]], left_on="pi_id", right_on="author", how="left"
    )

    data["publication_date"] = pd.to_datetime(data["publication_date"])

    return data


def preprocess_for_event_study(
    data: pd.DataFrame, applied_levels: pd.DataFrame
) -> pd.DataFrame:
    """
    Preprocesses the data for event study analysis.

    Args:
        data (pd.DataFrame): The input data to be preprocessed.
        applied_levels (pd.DataFrame): The applied_levels data used for merging.

    Returns:
        pd.DataFrame: The preprocessed data.

    """
    applied_levels["publication_date"] = pd.to_datetime(
        applied_levels["publication_date"]
    )
    applied_levels["parent_publication_date"] = pd.to_datetime(
        applied_levels["parent_publication_date"]
    )

    logger.info("Merging data with levels data")
    merged_data = pd.merge(
        data,
        applied_levels[["id", "parent_id", "parent_publication_date"]],
        on="id",
        how="inner",
    )

    logger.info("Sorting data by 'parent_id' and 'parent_publication_date'")
    merged_data.sort_values(
        by=["parent_id", "parent_publication_date"],
        ascending=[True, True],
        inplace=True,
    )

    logger.info("Dropping duplicates based on 'pi_id'")
    merged_data.drop_duplicates(subset="pi_id", keep="first", inplace=True)

    logger.info("Mapping 'parent_publication_date' to 'pi_id'")
    pi_id_date_dict = dict(
        zip(merged_data["pi_id"], merged_data["parent_publication_date"])
    )
    data["parent_publication_date"] = data["pi_id"].map(pi_id_date_dict)

    # data["parent_publication_date"] = "2022-01-01"

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
    Compute event study outputs based on the provided data and levels data.

    Args:
        data (pd.DataFrame): The input data.
        levels (pd.DataFrame): The levels data.

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

    data_pmid[["ca_publication_type", "ca_publication_date"]] = data_pmid[
        "cited_by_clin"
    ].apply(get_entrez_ptype_pmid)

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


def get_event_study_pc(data: pd.DataFrame, patent_data: pd.DataFrame):
    """
    Get event study data for protein concepts.

    Args:
        data (pd.DataFrame): The input data containing information about events.
        patent_data (pd.DataFrame): The patent data to be merged with the
        input data.

    Returns:
        tuple: A tuple containing three elements:
            - data (pd.DataFrame): The merged data.
            - final_data_counts (pd.DataFrame): The final data grouped by
            various columns and the count of occurrences.
            - final_data_citations (pd.DataFrame): The final data grouped
            by various columns and the sum of citations.
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
    final_data_citations = (
        final_data.groupby(["pi_id", "time", "seed", "protein_concept"])[
            "Cited by Patent Count"
        ]
        .sum()
        .reset_index(name="count")
    )

    return data, final_data_counts, final_data_citations


def get_applied_lab_staggered_outputs(
    data_loaders: AbstractDataset,
    mapping_df: AbstractDataset,
    sb_mapping_df: AbstractDataset,
    levels: pd.DataFrame,
    pdb_submissions: pd.DataFrame,
    strength_es: pd.DataFrame,
    patents_data: pd.DataFrame,
    clinical_citations: pd.DataFrame,
    grants_data: pd.DataFrame,
    mesh_terms: pd.DataFrame,
    institutional_data: pd.DataFrame,
):
    """
    Retrieves outputs from data loaders and performs data processing.

    Args:
        data_loaders (AbstractDataset): A dictionary-like object containing data loaders.
        mapping_df (AbstractDataset): A dataset containing mapping information.

    Returns:
        pd.DataFrame: The processed data.

    """
    # mesh terms extract tree_number first character
    mesh_terms["term_group"] = mesh_terms["tree_number"].apply(
        lambda x: str(x)[:1] if x is not None else None
    )

    # create dictionary of DUI to term_group
    mesh_terms_dict = mesh_terms.set_index("DUI")["term_group"].to_dict()

    # create dictionary of unique pi_id, strong using strength_es
    pi_strong = (
        strength_es.drop_duplicates(subset=["pi_id"])
        .set_index("pi_id")["strong"]
        .to_dict()
    )

    # change institutional_data columns to institution, unless column is institution
    institutional_data.drop_duplicates(subset="author", inplace=True)
    institutional_data.columns = [
        "institution_" + col if col != "institution" else col
        for col in institutional_data.columns
    ]

    mapping_df = mapping_df.loc[~mapping_df["author"].isin(sb_mapping_df["author"])]

    outputs = []
    agg_outputs = []
    collapsed_outputs = []
    for i, loader in enumerate(data_loaders.values()):
        logger.info("Loading data batch %d / %d", i + 1, len(data_loaders))
        data_batch = loader()

        # drop "display_name", "authorships", "ids" if they exist
        for col in ["display_name", "authorships", "ids"]:
            try:
                data_batch = data_batch.drop(columns=col)
            except KeyError:
                pass

        # drop rows with publication_date older than 2017-08-31
        data_batch = data_batch[data_batch["publication_date"] >= "2017-08-31"]

        # transform concepts to be a list of the ids (ie first item in sublists)
        data_batch["concepts"] = data_batch["concepts"].apply(
            lambda x: [y[0] for y in x] if x is not None else []
        )

        logger.info("Merging data with mapping information")
        mapping_df.drop_duplicates(subset="author", inplace=True)
        data_batch = data_batch.merge(
            mapping_df[["author", "seed"]],
            left_on="pi_id",
            right_on="author",
            how="left",
        )

        data_batch["publication_date"] = pd.to_datetime(data_batch["publication_date"])

        logger.info("Merging data with levels data")
        data_processed = _preprocess_for_staggered_design(data_batch, levels)

        logger.info("Merging data with pdb_submissions data")
        data_processed = _get_pdb_activity(data_processed, pdb_submissions)

        logger.info("Map strong links to labs")
        data_processed["strong"] = data_processed["pi_id"].map(pi_strong)
        data_processed.loc[data_processed["seed"] == "other", "strong"] = False

        logger.info("Map papers with protein_prediction links")
        data_processed["protein_concept"] = data_processed["concepts"].apply(
            lambda x: (
                any(element == "C18051474" or element == "C47701112" for element in x)
            )
        )

        logger.info("Map papers with experimental links")
        data_processed["experimental"] = data_processed["concepts"].apply(
            lambda x: (
                any(
                    element == "C185592680"
                    or element == "C55493867"
                    or element == "C12554922"
                    or element == "C46141821"
                    or element == "C121332964"
                    for element in x
                )
            )
        )

        logger.info("Create COVID 2020 share")
        data_processed = _collect_covid_references(data_processed)

        logger.info("Extracting fields")
        data_processed = _calculate_topic_share(data_processed)

        logger.info("Extracting mesh tags")
        data_processed = _calculate_mesh_balance(data_processed, mesh_terms_dict)

        logger.info("Merging data with patents data")
        data_processed = _get_patent_citations(data_processed, patents_data)

        logger.info("Merging data with clinical citations data")
        data_processed = _get_cc(data_processed, clinical_citations)

        logger.info("Merging data with grants data")
        data_processed = _get_awards(data_processed, grants_data)

        # drop columns
        data_processed.drop(
            columns=["concepts", "mesh_terms", "counts_by_year", "topics", "author"],
            inplace=True,
        )


        logger.info("Merging data with institutional information")
        data_processed = data_processed.merge(
            institutional_data,
            left_on="pi_id",
            right_on="institution_author",
            how="left",
        )

        # drop institution_author
        data_processed = data_processed.drop(columns=["institution_author"])

        quarterly, collapsed = _get_quarterly_aggregate_outputs(data_processed)
        # outputs.append(data_processed)
        agg_outputs.append(quarterly)
        collapsed_outputs.append(collapsed)

    # data = pd.concat(outputs)
    agg_data = pd.concat(agg_outputs)
    collapsed_data = pd.concat(collapsed_outputs)

    return agg_data, collapsed_data # TODO add bAack data


def _preprocess_for_staggered_design(
    data: pd.DataFrame, levels: pd.DataFrame
) -> pd.DataFrame:
    """
    Preprocesses the data for event study analysis.

    Args:
        data (pd.DataFrame): The input data to be preprocessed.
        levels (pd.DataFrame): The levels data used for merging.

    Returns:
        pd.DataFrame: The preprocessed data.

    """
    levels["publication_date"] = pd.to_datetime(levels["publication_date"])
    levels["parent_publication_date"] = pd.to_datetime(
        levels["parent_publication_date"]
    )

    logger.info("Merging data with levels data")
    merged_data = pd.merge(
        data,
        levels[["id", "parent_id", "parent_publication_date", "level"]],
        on="id",
        how="inner",
    )

    logger.info("Sorting data by 'parent_id' and 'parent_publication_date'")
    merged_data.sort_values(
        by=["parent_id", "parent_publication_date", "level"],
        ascending=[True, True, True],
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

    logger.info("Mapping 'parent_publication_date' and level to 'pi_id'")
    pi_id_date_dict = dict(
        zip(merged_data["pi_id"], merged_data["parent_publication_date"])
    )
    pi_id_level_dict = dict(zip(merged_data["pi_id"], merged_data["level"]))
    data["parent_publication_date"] = data["pi_id"].map(pi_id_date_dict)
    data["level"] = data["pi_id"].map(pi_id_level_dict)

    logger.info("Calculating yearly citations")
    data = _get_yearly_citations(data)

    data["before_af"] = data["publication_date"] < "2021-07-15"

    logger.info(
        "Computing time difference between 'publication_date' and 'parent_publication_date'"
    )
    data = _get_time(data)

    # Convert to integer
    final_data = data.copy()
    final_data.dropna(subset=["time"], inplace=True)
    final_data["time"] = final_data["time"].astype(int)

    return final_data
