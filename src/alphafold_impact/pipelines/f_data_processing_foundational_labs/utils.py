"""
Utility functions for data processing in the foundational labs pipeline.
"""

import logging
import pandas as pd
import numpy as np
from Bio import Entrez

Entrez.email = "david.ampudia@nesta.org.uk"

logger = logging.getLogger(__name__)

QUARTER_MAPPING = {  # HACK
    "2015Q3": -9,
    "2015Q4": -8,
    "2016Q1": -7,
    "2016Q2": -6,
    "2016Q3": -5,
    "2016Q4": -4,
    "2017Q1": -3,
    "2017Q2": -2,
    "2017Q3": -1,
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


def preprocess_for_staggered_design(
    data: pd.DataFrame, foundational_publications: pd.DataFrame
) -> pd.DataFrame:
    """
    Preprocesses the data for event study analysis.

    Args:
        data (pd.DataFrame): The input data to be preprocessed.
        foundational_publications (pd.DataFrame): The level0 data used for merging.

    Returns:
        pd.DataFrame: The preprocessed data.

    """
    foundational_publications["publication_date"] = pd.to_datetime(
        foundational_publications["publication_date"]
    )
    foundational_publications["parent_publication_date"] = pd.to_datetime(
        foundational_publications["parent_publication_date"]
    )

    logger.info("Merging data with level0 data")
    merged_data = pd.merge(
        data,
        foundational_publications[["id", "parent_id", "parent_publication_date"]],
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


def get_cum_sums(pi_data, publication_data):
    """
    Calculate the cumulative sums of counts for each PI in the given data.

    Args:
        pi_data (pd.DataFrame): The input data containing information about PIs.
        publication_data (pd.DataFrame): The publication data containing information
        about publications.

    Returns:
        pd.DataFrame: The data with cumulative sums of counts for each PI.
    """
    pi_relevant_pubs = publication_data.merge(
        pi_data[["id", "pi_id", "quarter"]], on="id", how="inner"
    )

    pi_relevant_pubs = (
        pi_relevant_pubs.groupby(["pi_id", "source", "quarter"])
        .size()
        .reset_index(name="counts")
    )
    pi_relevant_pubs = pi_relevant_pubs.sort_values(by=["pi_id", "source", "quarter"])
    pi_relevant_pubs["cumulative_counts"] = pi_relevant_pubs.groupby(
        ["pi_id", "source"]
    )["counts"].cumsum()

    pi_relevant_pubs = pi_relevant_pubs.pivot_table(
        index=["pi_id", "quarter"],
        columns="source",
        values="cumulative_counts",
        fill_value=0,
    ).reset_index()

    for col in ["af", "ct_ai", "ct_noai", "other"]:
        pi_relevant_pubs.rename(columns={col: "cum_" + col}, inplace=True)
        pi_relevant_pubs["cum_" + col] = pi_relevant_pubs["cum_" + col].astype(int)

    pi_data = pi_data.merge(pi_relevant_pubs, on=["pi_id", "quarter"], how="left")

    # sort by pi_id and quarter
    pi_data = pi_data.sort_values(by=["pi_id", "quarter"])

    # ffill the four columns
    pi_data[["cum_af", "cum_ct_ai", "cum_ct_noai", "cum_other"]] = (
        pi_data.groupby("pi_id")[["cum_af", "cum_ct_ai", "cum_ct_noai", "cum_other"]]
        .ffill()
        .fillna(0)
        .astype(int)
    )

    return pi_data


def get_strong_cum_sums(
    pi_data: pd.DataFrame, publication_data: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate the cumulative sums of strong counts for each PI in the given data.

    Args:
        pi_data (pd.DataFrame): The input data containing information about PIs.
        publication_data (pd.DataFrame): The publication data containing information

    Returns:
        pd.DataFrame: The data with cumulative sums of strong counts for each PI.
    """
    strong_publications = publication_data[
        publication_data["chain_label"].isin(["strong", "partial_strong"])
    ]

    pi_relevant_pubs = strong_publications.merge(
        pi_data[["id", "pi_id", "quarter"]], on="id", how="inner"
    )

    pi_relevant_pubs = (
        pi_relevant_pubs.groupby(["pi_id", "source", "quarter"])
        .size()
        .reset_index(name="counts")
    )

    pi_relevant_pubs = pi_relevant_pubs.sort_values(by=["pi_id", "source", "quarter"])

    pi_relevant_pubs["cumulative_counts"] = pi_relevant_pubs.groupby(
        ["pi_id", "source"]
    )["counts"].cumsum()

    pi_relevant_pubs = pi_relevant_pubs.pivot_table(
        index=["pi_id", "quarter"],
        columns="source",
        values="cumulative_counts",
        fill_value=0,
    ).reset_index()

    for col in ["af", "ct_ai", "ct_noai"]:
        pi_relevant_pubs.rename(columns={col: "strong_cumul_" + col}, inplace=True)
        pi_relevant_pubs["strong_cumul_" + col] = pi_relevant_pubs[
            "strong_cumul_" + col
        ].astype(int)

    pi_data = pi_data.merge(pi_relevant_pubs, on=["pi_id", "quarter"], how="left")

    # sort by pi_id and quarter
    pi_data = pi_data.sort_values(by=["pi_id", "quarter"])

    # ffill the four columns
    pi_data[
        [
            "strong_cumul_af",
            "strong_cumul_ct_ai",
            "strong_cumul_ct_noai",
        ]
    ] = (
        pi_data.groupby("pi_id")[
            [
                "strong_cumul_af",
                "strong_cumul_ct_ai",
                "strong_cumul_ct_noai",
            ]
        ]
        .ffill()
        .fillna(0)
        .astype(int)
    )

    return pi_data


def get_usage(data: pd.DataFrame) -> pd.DataFrame:
    """Get the usage of frontier and AF methods"""
    # get the authors
    author_data = (
        data.drop_duplicates(subset=["id"])
        .explode("authorships")
        .assign(
            author=lambda x: x["authorships"].apply(
                lambda y: y[0] if y is not None else None
            ),
            institution=lambda x: x["authorships"].apply(
                lambda y: y[1] if y is not None else None
            ),
            position=lambda x: x["authorships"].apply(
                lambda y: y[2] if y is not None else None
            ),
        )
        .dropna(subset=["author"])
        .drop_duplicates(subset=["id", "author"])
        .reset_index(drop=True)
        .drop(columns=["authorships"])
    )

    # drop A9999999999
    author_data = author_data[~(author_data["author"] == "A9999999999")]
    author_data["level"] = author_data["level"].astype(int)

    # create quarter column
    author_data["publication_date"] = pd.to_datetime(author_data["publication_date"])
    author_data["quarter"] = author_data["publication_date"].dt.to_period("Q")

    # sort by author, quarter, source
    author_data = author_data.sort_values(by=["author", "source", "quarter"])

    logger.info("Calculating cumulative counts")
    author_data = (
        author_data.groupby(["author", "source", "quarter"])
        .size()
        .reset_index(name="counts")
        .sort_values(by=["author", "source", "quarter"])
    )
    
    author_data["cumulative_counts"] = author_data.groupby(["author", "source"])["counts"].cumsum()

    logger.info("Pivoting table")
    author_data = author_data.pivot_table(
        index=["author", "quarter"],
        columns="source", 
        values="cumulative_counts",
        fill_value=0
    ).reset_index()

    # make source cols int
    for col in ["af", "ct_ai", "ct_noai", "other"]:
        author_data[col] = author_data[col].astype(int)

    return author_data



def get_intent(data: pd.DataFrame) -> pd.DataFrame:
    """
    Get the pi_id intent data from the given DataFrame.

    Args:
        data (pd.DataFrame): The input DataFrame containing the publications data.

    Returns:
        pd.DataFrame: A DataFrame containing the intent data.
    """

    # get the authors
    author_data = (
        data.drop_duplicates(subset=["id"])
        .explode("authorships")
        .assign(
            author=lambda x: x["authorships"].apply(
                lambda y: y[0] if y is not None else None
            ),
            institution=lambda x: x["authorships"].apply(
                lambda y: y[1] if y is not None else None
            ),
            position=lambda x: x["authorships"].apply(
                lambda y: y[2] if y is not None else None
            ),
        )
        .dropna(subset=["author"])
        .drop_duplicates(subset=["id", "author"])
        .reset_index(drop=True)
        .drop(columns=["authorships"])
    )

    # drop A9999999999
    author_data = author_data[~(author_data["author"] == "A9999999999")]

    # create quarter column
    author_data["publication_date"] = pd.to_datetime(author_data["publication_date"])
    author_data["quarter"] = author_data["publication_date"].dt.to_period("Q")

    # create grouped chain_labels
    labels = {
        "strong": ["strong", "partial_strong"],
        "weak": ["weak", "partial_weak"],
        "mixed": ["mixed"],
        "unknown": ["unknown", "no_data"],
    }

    # create author labels
    author_data["grouped_chain_label"] = author_data["chain_label"].apply(
        lambda x: next((k for k, v in labels.items() if x in v), "unknown")
    )

    author_labels = (
        author_data.groupby(["author", "quarter", "grouped_chain_label", "source"])
        .size()
        .reset_index(name="counts")
    )

    # pivot the author labels
    author_labels = author_labels.pivot_table(
        index=["author", "quarter"],
        columns=["source", "grouped_chain_label"],
        values="counts",
        fill_value=0,
    ).reset_index()

    # join the two-level columns into one
    author_labels.columns = [
        f"{col[0]}_{col[1]}" if isinstance(col, tuple) else col
        for col in author_labels.columns
    ]

    # rename author_ and quarter_ to author and quarter
    author_labels.rename(
        columns={"author_": "author", "quarter_": "quarter"}, inplace=True
    )

    # make all columns but author_ and quarter_ int
    author_labels = author_labels.astype(
        {
            col: "int"
            for col in author_labels.columns
            if col not in ["author", "quarter"]
        }
    )

    # create four columns that sum the other counts on row and substract the count of _unknown_
    author_labels["af_with_intent"] = (
        author_labels["af_strong"]
        + author_labels["af_weak"]
        + author_labels["af_mixed"]
    )
    author_labels["ct_ai_with_intent"] = (
        author_labels["ct_ai_strong"]
        + author_labels["ct_ai_weak"]
        + author_labels["ct_ai_mixed"]
    )
    author_labels["ct_noai_with_intent"] = (
        author_labels["ct_noai_strong"]
        + author_labels["ct_noai_weak"]
        + author_labels["ct_noai_mixed"]
    )
    author_labels["other_with_intent"] = (
        author_labels["other_strong"]
        + author_labels["other_weak"]
        + author_labels["other_mixed"]
    )

    count_cols = [
        col for col in author_labels.columns if col not in ["author", "quarter"]
    ]
    author_labels[count_cols] = author_labels.groupby("author")[count_cols].cumsum()

    return author_labels


def get_pdb_activity(data, pdb_submissions):
    """
    Calculate PDB activity metrics based on the given data and PDB submissions,
    adjusting the pdb_share to be the share of pdb submissions over the share of
    publications by each pi_id at each time.

    Args:
        data (pandas.DataFrame): The input data containing information about PIs,
          time, and publications count.
        pdb_submissions (pandas.DataFrame): The PDB submissions data containing
          information about PDB

    Returns:
        pandas.DataFrame: The updated data with additional PDB activity metrics.
    """
    # Convert resolution and R_free columns to numeric, coercing errors
    # Merge pdb_submissions with data on 'id'
    submissions = pdb_submissions.merge(data, on="id", how="inner")

    # Group by pi_id and time, then calculate metrics
    submissions_grouped = (
        submissions.groupby(["pi_id", "time"])
        .agg(
            num_uniprot_structures=pd.NamedAgg(
                column="num_uniprot_structures", aggfunc="sum"
            ),
            num_pdb_ids=pd.NamedAgg(column="num_pdb_ids", aggfunc="sum"),
            num_primary_submissions=pd.NamedAgg(
                column="num_primary_submissions", aggfunc="sum"
            ),
            score_mean=pd.NamedAgg(column="score_mean", aggfunc="mean"),
            complexity_sum=pd.NamedAgg(column="complexity_sum", aggfunc="sum"),
            complexity_mean=pd.NamedAgg(column="complexity_mean", aggfunc="mean"),
            organism_rarity_mean=pd.NamedAgg(
                column="organism_rarity_mean", aggfunc="mean"
            ),
            organism_rarity_max=pd.NamedAgg(
                column="organism_rarity_max", aggfunc="max"
            ),
            num_diseases=pd.NamedAgg(column="num_diseases", aggfunc="sum"),
            resolution_mean=pd.NamedAgg(column="resolution_mean", aggfunc="mean"),
            R_free_mean=pd.NamedAgg(column="R_free_mean", aggfunc="mean"),
            mean_tmscore=pd.NamedAgg(column="mean_tmscore", aggfunc="mean"),
            max_tmscore=pd.NamedAgg(column="max_tmscore", aggfunc="max"),
            normalised_mean_tmscore=pd.NamedAgg(
                column="normalised_mean_tmscore", aggfunc="mean"
            ),
            normalised_max_tmscore=pd.NamedAgg(
                column="normalised_max_tmscore", aggfunc="max"
            ),
        )
        .reset_index()
    )

    # Merge the grouped submissions back with the original data to include publication counts
    data_merged = data.merge(submissions_grouped, on=["pi_id", "time"], how="left")

    # get individual pdb_submission for each id
    data_merged = data_merged.merge(
        pdb_submissions[["id"]], on="id", how="left", indicator=True
    )
    data_merged["pdb_submission"] = data_merged["_merge"] == "both"
    data_merged.drop(columns=["_merge"], inplace=True)

    return data_merged


def get_patent_citations(data, patents_data):
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


def get_cc(data: pd.DataFrame, icite_data: pd.DataFrame) -> pd.DataFrame:
    """
    Merge clinical citations data with the given data and perform aggregation.

    Args:
        data (pandas.DataFrame): The input data.
        icite_data (pandas.DataFrame): The clinical citations data.

    Returns:
        pandas.DataFrame: The merged data with aggregated clinical citations information.
    """
    # preprocess icite data
    icite_data["pmid"] = icite_data["pmid"].astype(str)

    logger.info("Merging chains with iCite data")

    # merge on 'pmid'
    data_pmid = data.merge(icite_data[["pmid", "cited_by_clin"]], how="left", on="pmid")
    data_pmid = data_pmid.drop_duplicates(subset=["id"])

    logger.info("Exploding cited_by_clin")
    data_pmid = data_pmid[data_pmid["cited_by_clin"].astype(str) != "nan"]
    data_pmid = data_pmid[data_pmid["cited_by_clin"].notnull()]
    data_pmid["cited_by_clin"] = data_pmid["cited_by_clin"].apply(
        lambda x: x.split(" ")
    )
    data_pmid = data_pmid.explode("cited_by_clin")

    # reindex
    data_pmid.reset_index(inplace=True, drop=True)
    data_pmid.drop_duplicates(subset=["id", "cited_by_clin"], inplace=True)

    # merge clinical data on data
    data_merged = data[["id"]].merge(data_pmid, how="inner", on="id")

    # groupby id, get the count
    clinical_citations_grouped = (
        data_merged.groupby(["id"])
        .agg(
            ca_count=pd.NamedAgg(column="id", aggfunc="count"),
        )
        .reset_index()
    )

    # merge with data
    data = data.merge(clinical_citations_grouped, on="id", how="left")

    return data


def get_ai_use(concepts_list: list) -> bool:
    """
    Determines if any sublist in the provided list of concepts contains the AI
      concept identifier "C154945302".
    Args:
      concepts_list (list of lists): A list where each element is a sublist
        containing concept identifiers.
    Returns:
      bool: True if any sublist contains the AI concept identifier "C154945302",
        False otherwise.
    """

    # create a boolean if AI concept
    return any("C154945302" in sublist[0] for sublist in concepts_list)


def get_awards(data, grants_data):
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


def calculate_field_share(data):
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


def calculate_primary_field(data):
    """
    Processes the input DataFrame to calculate the primary field for each author and time,
    and returns the modified DataFrame with additional columns representing the count of
    each primary field.
    Args:
        data (pd.DataFrame): Input DataFrame containing columns 'author', 'time', and 'topics'.
            'topics' should be a list of lists where each sublist contains
            topic information.
    Returns:
        pd.DataFrame: Modified DataFrame with additional columns for each primary field count,
            prefixed with 'primary_field_'.
    """

    df = data.copy()

    # groupby and calculate count
    df = (
        df.groupby(["author", "time", "primary_field"]).size().reset_index(name="count")
    )

    # pivot the DataFrame to get one column for each primary field
    df = df.pivot(index=["author", "time"], columns="primary_field", values="count")

    # reset index
    df.reset_index(inplace=True)

    # fill NaN values with 0
    df = df.fillna(0)

    # change column names to be camel case, and prefix with "field_"
    df.columns = [
        "primary_field_" + col if col != "author" and col != "time" else col
        for col in df.columns
    ]

    # merge with data on author and time
    data = data.merge(df, on=["author", "time"], how="left")

    return data


def calculate_mesh_balance(
    data: pd.DataFrame, mesh_terms: pd.DataFrame
) -> pd.DataFrame:
    """
    Calculate the balance of mesh terms for each author over time.

        Args:
        data (pd.DataFrame): A DataFrame containing author data with columns 'author',
          'time', and 'mesh_terms'.
        mesh_terms (pd.DataFrame): A DataFrame containing mesh terms and their
          corresponding groups.
    Returns:
        pd.DataFrame: A DataFrame with the original data and additional columns for
          each mesh term's share, prefixed with 'mesh_'.
    """

    df = data.copy()

    # transform mesh_terms by mapping to term_group
    df["mesh_terms"] = df["mesh_terms"].apply(
        lambda x: [y[0] for y in x] if x is not None else []
    )
    df["mesh_terms"] = df["mesh_terms"].apply(
        lambda x: [mesh_terms.get(y, "") for y in x] if x is not None else []
    )

    # explode the mesh_terms column into multiple rows
    df = df.explode("mesh_terms")

    # group by author and mesh_terms and calculate the count
    df = df.groupby(["author", "quarter", "mesh_terms"]).size().reset_index(name="count")

    # calculate the total count for each author
    total_count = df.groupby("author")["count"].sum()

    # calculate the share of each mesh term
    df["share"] = df.apply(
        lambda row: row["count"] / total_count[row["author"]], axis=1
    )

    # pivot the DataFrame to get one column for each mesh term
    df = df.pivot(index=["author", "quarter"], columns="mesh_terms", values="share")

    # reset index
    df.reset_index(inplace=True)

    # change column names to be camel case, and prefix with "mesh_"
    df.columns = [
        "mesh_" + col if col != "author" and col != "quarter" else col
        for col in df.columns
    ]

    # fill NaN values with 0
    df = df.fillna(0)

    # try to drop "mesh_" column
    if "mesh_" in df.columns:
        df = df.drop(columns=["mesh_"])

    # merge with data on author and time
    data = data.merge(df, on=["author", "quarter"], how="left")

    data.drop(columns=["mesh_terms"], inplace=True)

    # fill NaN for columns that are mesh_
    mesh_cols = [col for col in data.columns if col.startswith("mesh_")]
    data[mesh_cols] = data[mesh_cols].fillna(0)
    
    return data


def collect_covid_references(data: pd.DataFrame) -> pd.DataFrame:
    """
    Collects and calculates the share of COVID-19 related concepts for each principal
     investigator (pi_id) in the provided DataFrame for the months of September to
     December 2020.

    Args:
      data (pd.DataFrame): A DataFrame containing at least the following columns:
    Returns:
      pd.DataFrame: The original DataFrame with an additional column "covid_share_2020"
    """

    data_2020 = data[
        (data["publication_date"] >= "2020-01-01")
        & (data["publication_date"] <= "2020-12-31")
    ]

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


def get_quarterly_aggregate_outputs(data):
    """
    Aggregate the data by 'pi_id' and 'time' (quarter) and calculate various metrics.

    Args:
        data (pandas.DataFrame): The input data containing the columns 'pi_id', 'time',
        and other metrics.

    Returns:
        pandas.DataFrame: The aggregated data with calculated metrics.

    """

    def safe_mode(series):
        mode = series.mode()
        return mode.iloc[0] if not mode.empty else np.nan

    # Define aggregation functions for each column
    agg_funcs = {
        "cited_by_count": "sum",
        "seed": "first",
        "cit_0": "sum",
        "cit_1": "sum",
        "parent_time": "first",
        "num_uniprot_structures": safe_mode,
        "num_pdb_ids": safe_mode,
        "num_primary_submissions": safe_mode,
        "score_mean": safe_mode,
        "complexity_sum": safe_mode,
        "complexity_mean": safe_mode,
        "organism_rarity_mean": safe_mode,
        "organism_rarity_max": safe_mode,
        "num_diseases": safe_mode,
        "resolution_mean": safe_mode,
        "R_free_mean": safe_mode,
        "mean_tmscore": safe_mode,
        "max_tmscore": safe_mode,
        "normalised_mean_tmscore": safe_mode,
        "normalised_max_tmscore": safe_mode,
        "pdb_submission": "sum",
        "intent": "first",
        "ai_concept": lambda x: x.sum(),
        "protein_concept": lambda x: x.sum(),
        "experimental": lambda x: x.sum(),
        "covid_share_2020": "first",
        "patent_count": "sum",
        "patent_citation": "sum",
        "CPCs": lambda x: ";;".join(
            filter(lambda i: i != "nan" and i != "None", map(str, x))
        ),
        "ca_count": "sum",
        "fwci": "mean",
        "citation_normalized_percentile_value": "mean",
        "citation_normalized_percentile_is_in_top_1_percent": lambda x: int(x.sum()),
        "citation_normalized_percentile_is_in_top_10_percent": lambda x: int(x.sum()),
        "cum_af": safe_mode,
        "cum_ct_ai": safe_mode,
        "cum_ct_noai": safe_mode,
        "cum_other": safe_mode,
        "strong_cumul_af": safe_mode,
        "strong_cumul_ct_ai": safe_mode,
        "strong_cumul_ct_noai": safe_mode,
        "primary_field": safe_mode,
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

    # remove ";;None" in CPCs
    data_grouped["CPCs"] = data_grouped["CPCs"].str.replace(";;None", "")
    data_grouped["CPCs"] = data_grouped["CPCs"].str.replace("None;;", "")

    # adding counts
    data_grouped["num_publications"] = data.groupby(["pi_id", "time"]).size().values

    return data_grouped


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

    # Replace the quarters with the corresponding numbers
    data["time"] = data["publication_date"].map(QUARTER_MAPPING)
    data["parent_time"] = data["parent_publication_date"].map(QUARTER_MAPPING)

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

    t_columns.columns = ["cit_" + str(int(i)) for i in t_columns.columns]

    # fill remaining NaN values in 't0' with 0
    t_columns = t_columns.fillna(0)
    t_columns.reset_index(inplace=True)

    # merge the 't' DataFrame back to the original DataFrame
    data = data.merge(t_columns[["id", "cit_0", "cit_1"]], on="id", how="left")

    # fill ct0 and ct1 with nan
    data["cit_0"] = data["cit_0"].fillna(0)
    data["cit_1"] = data["cit_1"].fillna(0)

    return data


def process_institutional_data(institutional_data: pd.DataFrame) -> pd.DataFrame:
    """Process the institutional data by dropping duplicates and renaming columns."""
    # change institutional_data columns to institution, unless column is institution
    institutional_data.drop_duplicates(subset="institution", inplace=True)
    institutional_data.columns = [
        "institution_" + col if col != "institution" else col
        for col in institutional_data.columns
    ]

    return institutional_data


def process_pdb_data(pdb_submissions: pd.DataFrame) -> pd.DataFrame:
    """Process the PDB submissions data by exploding the authorships."""
    return (
        pdb_submissions.explode("authorships")
        .assign(
            authorships=lambda x: x["authorships"].apply(
                lambda y: y[0] if isinstance(y, np.ndarray) and len(y) > 0 else None
            )
        )
        .rename(columns={"authorships": "author"})
    )
