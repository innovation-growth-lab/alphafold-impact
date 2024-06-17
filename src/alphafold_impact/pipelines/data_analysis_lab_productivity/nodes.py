"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
from kedro.io import AbstractDataset
from joblib import Parallel, delayed
from Bio import Entrez
from ..data_analysis_descriptive_translational.nodes import (  # pylint: disable=relative-beyond-top-level
    get_entrez_ptype_pmid,
)


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

        # drop "display_name", "authorships", "ids" if they exist
        for col in ["display_name", "authorships", "ids"]:
            try:
                data_batch = data_batch.drop(columns=col)
            except KeyError:
                pass

        # drop rows with publication_date older than 2017-08-31
        data_batch = data_batch[data_batch["publication_date"] >= "2017-08-31"]

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
        patent_data (pd.DataFrame): The patent data to be merged with the input data.

    Returns:
        tuple: A tuple containing three elements:
            - data (pd.DataFrame): The merged data.
            - final_data_counts (pd.DataFrame): The final data grouped by various columns
                and the count of occurrences.
            - final_data_citations (pd.DataFrame): The final data grouped by various columns
                and the sum of citations.
    """

    # merge patent_data on data
    data = data.merge(
        patent_data, how="inner", left_on="doi", right_on="NPL Resolved External ID(s)"
    )

    final_data = data.copy()
    final_data.dropna(subset=["time"], inplace=True)
    final_data["time"] = final_data["time"].astype(int)

    final_data_counts = (
        final_data.groupby(["pi_id", "time", "seed"]).size().reset_index(name="count")
    )
    final_data_citations = (
        final_data.groupby(["pi_id", "time", "seed"])["Cited by Patent Count"]
        .sum()
        .reset_index(name="count")
    )

    return data, final_data_counts, final_data_citations


def get_sb_lab_staggered_outputs(
    data_loaders: AbstractDataset,
    mapping_df: AbstractDataset,
    level0: pd.DataFrame,
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

    # extract author id (first item in each sublist in authorships)
    level0_authors = level0.copy()
    level0_authors["authorships"] = level0["authorships"].apply(
        lambda x: [y[0] for y in x] if x is not None else []
    )
    level0_authors = level0_authors.explode("authorships")
    level0_authors.drop_duplicates(subset=["id", "authorships"], inplace=True)

    # create a dictionary for each author, with "af", "ct", "other" counts (the values in column source)
    author_counts = (
        level0_authors.groupby(["authorships", "source"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    author_counts.columns = [
        "ext_" + col if col != "authorships" else col for col in author_counts.columns
    ]

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

        logger.info("Merging data with level0 data")
        data_processed = _preprocess_for_staggered_design(data_batch, level0)

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

        logger.info("Collect COVID 2020 references")
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

        # add extensive margin
        data_processed = data_processed.merge(
            author_counts, left_on="pi_id", right_on="authorships", how="left"
        )

        # drop topics, mesh_terms, concepts
        data_processed = data_processed.drop(
            columns=["author", "topics", "mesh_terms", "concepts", "counts_by_year"]
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
        agg_outputs.append(quarterly)
        collapsed_outputs.append(collapsed)

    agg_data = pd.concat(agg_outputs)
    collapsed_data = pd.concat(collapsed_outputs)

    return agg_data, collapsed_data


def _preprocess_for_staggered_design(
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


def _get_time(data):
    """
    Converts publication dates to quarters and maps them to corresponding numbers.

    Args:
        data (pandas.DataFrame): The input data containing publication dates.

    Returns:
        pandas.DataFrame: The input data with publication dates converted to quarters
        and mapped to numbers.
    """

    data["publication_date"] = data["publication_date"].dt.to_period("Q").astype(str)
    data["parent_publication_date"] = (
        data["parent_publication_date"].dt.to_period("Q").astype(str)
    )

    # if seed isother, set parent_publication_date to NaT
    data.loc[data["seed"] == "other", "parent_publication_date"] = pd.NaT

    # Create a mapping dictionary
    quarter_mapping = {
        "2017Q4": 0,
        "2018Q1": 1,
        "2018Q2": 2,
        "2018Q3": 3,
        "2018Q4": 4,
        "2019Q1": 5,
        "2019Q2": 6,
        "2019Q3": 7,
        "2019Q4": 8,
        "2020Q1": 9,
        "2020Q2": 10,
        "2020Q3": 11,
        "2020Q4": 12,
        "2021Q1": 13,
        "2021Q2": 14,
        "2021Q3": 15,
        "2021Q4": 16,
        "2022Q1": 17,
        "2022Q2": 18,
        "2022Q3": 19,
        "2022Q4": 20,
        "2023Q1": 21,
        "2023Q2": 22,
        "2023Q3": 23,
        "2023Q4": 24,
        "2024Q1": 25,
    }

    # Replace the quarters with the corresponding numbers
    data["time"] = data["publication_date"].map(quarter_mapping)
    data["parent_time"] = data["parent_publication_date"].map(quarter_mapping)

    return data


def _get_yearly_citations(data):
    """
    Calculate yearly citations for each publication in the given data.

    Args:
        data (pd.DataFrame): The input data containing publication information.

    Returns:
        pd.DataFrame: The input data with additional columns for yearly citations.

    """

    data_exploded = data.explode("counts_by_year")
    data_exploded.reset_index(drop=True, inplace=True)

    data_exploded[["cited_by_count", "year"]] = pd.json_normalize(
        data_exploded["counts_by_year"]
    )

    data_exploded["publication_date"] = pd.to_datetime(
        data_exploded["publication_date"]
    ).dt.year

    # calculate the 't' value
    data_exploded["ct"] = data_exploded["year"] - data_exploded["publication_date"]

    data_exploded.dropna(subset=["ct"], inplace=True)
    data_exploded["ct"] = data_exploded["ct"].astype(int)

    t_columns = data_exploded.pivot_table(
        index="id", columns="ct", values="cited_by_count", aggfunc="sum"
    )

    # remove negative columns
    t_columns = t_columns.loc[:, (t_columns.columns >= 0)]

    t_columns.columns = ["ct" + str(int(i)) for i in t_columns.columns]

    # fill remaining NaN values in 't0' with 0
    t_columns = t_columns.fillna(0)
    t_columns.reset_index(inplace=True)

    # merge the 't' DataFrame back to the original DataFrame
    data = data.merge(t_columns[["id", "ct0", "ct1", "ct2"]], on="id", how="left")

    # fill ct0 and ct1 with nan
    data["ct0"] = data["ct0"].fillna(0)
    data["ct1"] = data["ct1"].fillna(0)

    return data


def _get_pdb_activity(data, pdb_submissions):
    """
    Calculate PDB activity metrics based on the given data and PDB submissions.

    Args:
        data (pandas.DataFrame): The input data containing information about PIs
          and time.
        pdb_submissions (pandas.DataFrame): The PDB submissions data containing
          information about PDB IDs, resolution, and R_free.

    Returns:
        pandas.DataFrame: The updated data with additional PDB activity metrics.

    """
    pdb_submissions["resolution"] = pd.to_numeric(
        pdb_submissions["resolution"], errors="coerce"
    )
    pdb_submissions["R_free"] = pd.to_numeric(
        pdb_submissions["R_free"], errors="coerce"
    )

    submissions = pdb_submissions[["id", "resolution", "R_free"]].merge(
        data, on="id", how="inner"
    )

    total_count = submissions["id"].count()
    submissions_grouped = (
        submissions.groupby(["pi_id", "time"])
        .agg(
            pdb_share=pd.NamedAgg(column="id", aggfunc=lambda x: len(x) / total_count),
            resolution=pd.NamedAgg(column="resolution", aggfunc="mean"),
            R_free=pd.NamedAgg(column="R_free", aggfunc="mean"),
        )
        .reset_index()
    )

    # merge with data
    data = data.merge(submissions_grouped, on=["pi_id", "time"], how="left")

    return data


def _get_patent_citations(data, patents_data):
    """
    Get patent citations for the given data.

    Args:
        data (pandas.DataFrame): The input data.
        patents_data (pandas.DataFrame): The patent data.

    Returns:
        pandas.DataFrame: The data with patent citations merged and aggregated.
    """

    # merge patent_data on data
    patent_matches = data.merge(
        patents_data, how="inner", left_on="doi", right_on="NPL Resolved External ID(s)"
    )

    # groupby id, get the count and the joined of CPCs
    patent_matches_grouped = (
        patent_matches.groupby(["id", "time"])
        .agg(
            patent_count=pd.NamedAgg(column="id", aggfunc="count"),
            CPCs=pd.NamedAgg(
                column="CPC Classifications",
                aggfunc=lambda x: ";;".join(map(str, x)),
            ),
            patent_citation=pd.NamedAgg(column="Cited by Patent Count", aggfunc="sum"),
        )
        .reset_index()
    )

    # merge with data
    data = data.merge(patent_matches_grouped, on=["id", "time"], how="left")

    # fill count with 0
    data["patent_count"] = data["patent_count"].fillna(0)

    return data


def _get_cc(data, clinical_citations):
    """
    Merge clinical citations data with the given data and perform aggregation.

    Args:
        data (pandas.DataFrame): The input data.
        clinical_citations (pandas.DataFrame): The clinical citations data.

    Returns:
        pandas.DataFrame: The merged data with aggregated clinical citations information.
    """

    # merge clinical data on data
    data_merged = data[["id"]].merge(clinical_citations, how="inner", on="id")
    data_merged["ca_publication_date"] = pd.to_datetime(
        data_merged["ca_publication_date"]
    )

    # compute the time difference between publication_date and ca_publication_date
    data_merged["tcc"] = (
        data_merged["ca_publication_date"] - data_merged["publication_date"]
    ) / pd.Timedelta(days=90)

    # groupby id, get the count and joined ca_publication_type. Get the average ca_publciation_date
    clinical_citations_grouped = (
        data_merged.groupby(["id"])
        .agg(
            ca_count=pd.NamedAgg(column="id", aggfunc="count"),
            ca_publication_type=pd.NamedAgg(
                column="ca_publication_type", aggfunc=lambda x: ";;".join(map(str, x))
            ),
            ca_publication_date=pd.NamedAgg(
                column="ca_publication_date", aggfunc="mean"
            ),
            tcc=pd.NamedAgg(column="tcc", aggfunc="mean"),
        )
        .reset_index()
    )

    # fix format of ca_publication_date
    clinical_citations_grouped["ca_publication_date"] = (
        pd.to_datetime(clinical_citations_grouped["ca_publication_date"])
        .dt.to_period("Q")
        .astype(str)
    )

    # merge with data
    data = data.merge(clinical_citations_grouped, on="id", how="left")

    return data


def _get_awards(data, grants_data):
    """
    Retrieves awards information from the given data and grants_data.

    Args:
        data (pandas.DataFrame): The input data containing information about awards.
        grants_data (pandas.DataFrame): The grants data containing additional information
        about awards.

    Returns:
        pandas.DataFrame: The merged data with awards information.

    """
    # merge grants_data on data
    grant_matches = data.merge(grants_data, how="inner", on="id")

    # drop "id", "funder_display_name" duplicates
    grant_matches = grant_matches.drop_duplicates(subset=["id", "funder_display_name"])

    # groupby id, get the count and the joined of CPCs
    grant_matches_grouped = (
        grant_matches.groupby(["id"])
        .agg(
            grant_count=pd.NamedAgg(column="id", aggfunc="count"),
            grant_agency=pd.NamedAgg(
                column="funder_display_name", aggfunc=lambda x: ";;".join(map(str, x))
            ),
        )
        .reset_index()
    )

    # merge with data
    data = data.merge(grant_matches_grouped, on=["id"], how="left")

    # fill count with 0
    data["grant_count"] = data["grant_count"].fillna(0)

    return data


def _calculate_topic_share(data):
    """
    Calculate the share of each subfield for each author in a DataFrame.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data.

    Returns:
        pandas.DataFrame: The resulting DataFrame with the share of each subfield for each author.
    """
    df = data.copy()
    df["fields"] = df["topics"].apply(
        lambda x: [y[5] for y in x] if x is not None and len(x) > 0 else []
    )
    # explode the subfields column into multiple rows
    df = df.explode("fields")
    # group by author and subfields and calculate the count
    df = df.groupby(["author", "time", "fields"]).size().reset_index(name="count")
    # calculate the total count for each author
    total_count = df.groupby("author")["count"].sum()
    # calculate the share of each subfield
    df["share"] = df.apply(
        lambda row: row["count"] / total_count[row["author"]], axis=1
    )
    # pivot the DataFrame to get one column for each subfield
    df = df.pivot(index=["author", "time"], columns="fields", values="share")

    # reset index
    df.reset_index(inplace=True)

    # fill NaN values with 0
    df = df.fillna(0)

    # change column names to be camel case, and prefix with "field_"
    df.columns = [
        (
            "field_" + col.lower().replace(" ", "_")
            if col != "author" and col != "time"
            else col
        )
        for col in df.columns
    ]

    # merge with data on author and time
    data = data.merge(df, on=["author", "time"], how="left")

    return data


def _calculate_mesh_balance(data, mesh_terms_dict):
    df = data.copy()

    # transform mesh_terms by mapping to term_group
    df["mesh_terms"] = df["mesh_terms"].apply(
        lambda x: [y[0] for y in x] if x is not None else []
    )
    df["mesh_terms"] = df["mesh_terms"].apply(
        lambda x: [mesh_terms_dict.get(y, "") for y in x] if x is not None else []
    )
    # iterate over mesh_terms column values, extracting the first character of each element in the nested list

    # explode the mesh_terms column into multiple rows
    df = df.explode("mesh_terms")

    # group by author and mesh_terms and calculate the count
    df = df.groupby(["author", "time", "mesh_terms"]).size().reset_index(name="count")

    # calculate the total count for each author
    total_count = df.groupby("author")["count"].sum()

    # calculate the share of each mesh term
    df["share"] = df.apply(
        lambda row: row["count"] / total_count[row["author"]], axis=1
    )

    # pivot the DataFrame to get one column for each mesh term
    df = df.pivot(index=["author", "time"], columns="mesh_terms", values="share")

    # reset index
    df.reset_index(inplace=True)

    # change column names to be camel case, and prefix with "mesh_"
    df.columns = [
        "mesh_" + col if col != "author" and col != "time" else col
        for col in df.columns
    ]

    # fill NaN values with 0
    df = df.fillna(0)

    # merge with data on author and time
    data = data.merge(df, on=["author", "time"], how="left")

    return data


def _collect_covid_references(data):
    data_2020 = data[data["time"].isin([9, 10, 11, 12])]

    data_2020["covid_2020"] = data_2020["concepts"].apply(
        lambda x: (any(element == "C524204448" for element in x))
    )

    # Group the filtered dataframe by pi_id and calculate the share of COVID concepts
    covid_share_2020 = (
        data_2020.groupby("pi_id")["covid_2020"].sum()
        / data_2020.groupby("pi_id")["covid_2020"].count()
    )

    # Convert the Series to a DataFrame
    covid_share_2020 = covid_share_2020.reset_index()

    # Rename the column
    covid_share_2020.columns = ["pi_id", "covid_share_2020"]

    # Create a dictionary from the DataFrame
    covid_share_dict = dict(
        zip(covid_share_2020["pi_id"], covid_share_2020["covid_share_2020"])
    )

    data["covid_share_2020"] = data["pi_id"].map(covid_share_dict)

    return data


def _get_quarterly_aggregate_outputs(data):
    """
    Aggregate the data by 'pi_id' and 'time' (quarter) and calculate various metrics.

    Args:
        data (pandas.DataFrame): The input data containing the columns 'pi_id', 'time',
        and other metrics.

    Returns:
        pandas.DataFrame: The aggregated data with calculated metrics.

    """

    # Define aggregation functions for each column
    agg_funcs = {
        "cited_by_count": "sum",
        "seed": "first",
        "ct0": "sum",
        "ct1": "sum",
        "tcc": "mean",
        "parent_time": "first",
        "pdb_share": "first",
        "resolution": "mean",
        "R_free": "mean",
        "strong": "first",
        "protein_concept": lambda x: x.sum(),
        "experimental": lambda x: x.sum(),
        "covid_share_2020": "first",
        "patent_count": "sum",
        "patent_citation": "sum",
        "CPCs": lambda x: ";;".join(
            filter(lambda i: i != "nan" and i != "None", map(str, x))
        ),
        "ca_count": "sum",
        "ca_publication_type": lambda x: ";;".join(
            filter(lambda i: i != "nan" and i != "None", map(str, x))
        ),
        "ca_publication_date": "first",
        "grant_count": "sum",
        "grant_agency": lambda x: ";;".join(
            filter(lambda i: i != "nan" and i != "None", map(str, x))
        ),
    }

    # Assume `data` is your original data
    for column in data.columns:
        if (
            column.startswith("field_")
            or column.startswith("mesh_")
            or column.startswith("institution_")
            or column.startswith("ext_")
        ):
            agg_funcs[column] = "first"

    # Group by 'pi_id' and 'time' (quarter) and aggregate
    data_grouped = data.groupby(["pi_id", "time"]).agg(agg_funcs).reset_index()

    # remove ";;None" in CPCs, ca_publication_type, grant_agency
    for col in ["CPCs", "ca_publication_type", "grant_agency"]:
        data_grouped[col] = data_grouped[col].str.replace(";;None", "")
        data_grouped[col] = data_grouped[col].str.replace("None;;", "")

    # adding counts
    data_grouped["num_publications"] = data.groupby(["pi_id", "time"]).size().values

    # create collapsed data
    data_collapsed = data.groupby(["pi_id", "before_af"]).agg(agg_funcs).reset_index()

    # remove ";;None" in CPCs, ca_publication_type, grant_agency
    for col in ["CPCs", "ca_publication_type", "grant_agency"]:
        data_collapsed[col] = data_collapsed[col].str.replace(";;None", "")
        data_collapsed[col] = data_collapsed[col].str.replace("None;;", "")

    # adding counts
    data_collapsed["num_publications"] = (
        data.groupby(["pi_id", "before_af"]).size().values
    )

    return data_grouped, data_collapsed
