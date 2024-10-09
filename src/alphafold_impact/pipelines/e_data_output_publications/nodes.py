"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd

logger = logging.getLogger(__name__)


def create_publications_data(
    data: pd.DataFrame,
    source: str,
    mesh_terms: pd.DataFrame,
    patents_data,
    pdb_submissions,
    icite_data,
):
    """
    Processes and enriches publication data with various metrics and annotations.
    Args:
        data (pd.DataFrame): The primary data containing publication information.
        source (str): The source of the data.
        mesh_terms (pd.DataFrame): DataFrame containing MeSH terms and their tree numbers.
        patents_data: Data related to patent citations.
        pdb_submissions: Data related to PDB submissions.
        icite_data: Data from the iCite database.
    Returns:
        pd.DataFrame: The enriched publication data with additional metrics and annotations.
    """
    if "level" in data.columns:
        data = data[data["level"] != 3]
    try:
        data = _sort_drop(data)
    except Exception as e:  # pylint: disable=broad-except
        logger.info(e)
        data.drop_duplicates(subset="id", inplace=True)

    mesh_terms["term_group"] = mesh_terms["tree_number"].apply(
        lambda x: str(x)[:1] if x is not None else None
    )
    mesh_terms_dict = mesh_terms.set_index("DUI")["term_group"].to_dict()

    data["source"] = source

    # get fields data
    data = _get_field_distr(data)

    # get mesh terms
    data = _calculate_mesh_balance(data, mesh_terms_dict)

    # fill columns that begin with "mesh_" with 0
    data[[col for col in data.columns if col.startswith("mesh_")]] = data[
        [col for col in data.columns if col.startswith("mesh_")]
    ].fillna(0)

    # get patent citations
    data = _get_patent_citations(data, patents_data)

    # add pdb activity metrics
    data = _get_pdb_activity(data, pdb_submissions)

    # drop mesh_terms column
    data.drop(columns="mesh_terms", inplace=True)

    # icite data
    icite_outputs = _create_cc_counts(
        data[["parent_id", "id", "level", "source", "pmid", "doi"]], icite_data
    )

    data = data.merge(icite_outputs[["id", "count"]], on="id", how="left")

    return data


def _sort_drop(data):
    sort_order = {
        "strong": 0,
        "partial_strong": 1,
        "mixed": 2,
        "weak": 3,
        "partial_weak": 4,
        "unknown": 5,
        "no_data": 6,
    }
    data["sort_order"] = data["chain_label"].map(sort_order)

    data = data.sort_values("sort_order").groupby(["level", "id"]).first().reset_index()

    data = data[
        ~data["id"].isin(["W3177828909", "W3211795435", "W3202105508", "W3202105508"])
    ]

    data = data[
        ~(
            (data["parent_id"].isin(["W3177828909", "W3211795435", "W3202105508"]))
            & (data["level"] != 0)
        )
    ]

    data.sort_values("sort_order").drop_duplicates(subset="id", inplace=True)

    # drop the sort_order column, reindex
    data.drop(columns="sort_order", inplace=True)
    data.reset_index(drop=True, inplace=True)
    return data


def _create_cc_counts(data, icite_data):
    """
    Creates a count column in the given data DataFrame based on the 'cited_by_clin' column
        in the icite_data DataFrame.

    Args:
        data (pandas.DataFrame): The input data DataFrame.
        icite_data (pandas.DataFrame): The icite_data DataFrame containing the
            'cited_by_clin' column.

    Returns:
        pandas.DataFrame: The updated data DataFrame with a new 'count' column.
    """
    # change pmid, doi to str
    icite_data["pmid"] = icite_data["pmid"].astype(str)
    icite_data["doi"] = icite_data["doi"].astype(str)

    # merge on 'pmid'
    data_pmid = data.merge(icite_data[["pmid", "cited_by_clin"]], how="left", on="pmid")
    data_pmid = data_pmid.drop_duplicates(subset=["parent_id", "id", "level", "source"])

    # merge on 'doi'
    data_doi = data.merge(icite_data[["doi", "cited_by_clin"]], how="left", on="doi")
    data_doi = data_doi.drop_duplicates(subset=["parent_id", "id", "level", "source"])

    combined_data = pd.concat([data_pmid, data_doi]).drop_duplicates(
        subset=["parent_id", "id", "level", "source"]
    )

    logger.info("Exploding cited_by_clin")
    combined_data = combined_data[combined_data["cited_by_clin"].astype(str) != "nan"]

    # drop if cited_by_clin is None
    combined_data = combined_data[combined_data["cited_by_clin"].notnull()]

    combined_data["cited_by_clin"] = combined_data["cited_by_clin"].apply(
        lambda x: x.split(" ")
    )

    # create count column
    combined_data["count"] = combined_data["cited_by_clin"].apply(len)

    # merge back with data, fill with 0 for missing count
    data = data.merge(combined_data[["id", "count"]], on="id", how="left")
    data["count"] = data["count"].fillna(0)

    return data


def merge_individual_data(
    data_af: pd.DataFrame, data_ct: pd.DataFrame, data_other: pd.DataFrame
) -> pd.DataFrame:
    """
    Merges individual dataframes by removing overlapping IDs and concatenating them.

    Args:
        data_af (pd.DataFrame): DataFrame containing AlphaFold data.
        data_ct (pd.DataFrame): DataFrame containing CT data.
        data_other (pd.DataFrame): DataFrame containing other data.

    Returns:
        pd.DataFrame: Merged DataFrame with cleaned and concatenated data.
    """

    # drop id in data_ct if it is in data_af
    data_ct = data_ct[~data_ct["id"].isin(data_af["id"])]

    # drop id in data_other if it is in data_af or data_ct
    data_other = data_other[~data_other["id"].isin(data_af["id"])]

    # concatenate, ignore index
    data = pd.concat([data_af, data_ct, data_other], ignore_index=True)

    # drop column '0' if it exists
    if "0" in data.columns:
        data.drop(columns="0", inplace=True)

    # change column level to str, drop -1
    data["level"] = data["level"].astype(str)
    data = data[data["level"] != "-1"]

    return data


def _get_field_distr(data: pd.DataFrame) -> pd.DataFrame:
    """
    Processes the input DataFrame to calculate the distribution of fields for each author.
    Args:
        data (pd.DataFrame): Input DataFrame containing at least the columns 'id' and 'topics'.
            'topics' should be a list of lists where each sublist contains elements
            with the 6th element being the field of interest.
    Returns:
        pd.DataFrame: A DataFrame with the original data merged with the calculated field
            distribution. The resulting DataFrame will have additional columns prefixed
            with 'field_' representing the share of each field for each author.
    """

    df = data.copy()
    df["fields"] = df["topics"].apply(
        lambda x: [y[5] for y in x] if x is not None and len(x) > 0 else []
    )
    # explode the subfields column into multiple rows
    df = df.explode("fields")
    # group by author and subfields and calculate the count
    df = df.groupby(["id", "fields"]).size().reset_index(name="count")

    # calculate the share of each subfield
    df["share"] = df.apply(lambda row: row["count"] / 1, axis=1)
    # pivot the DataFrame to get one column for each subfield
    df = df.pivot(index=["id"], columns="fields", values="share")

    # reset index
    df.reset_index(inplace=True)

    # fill NaN values with 0
    df = df.fillna(0)

    # change column names to be camel case, and prefix with "field_"
    df.columns = [
        ("field_" + col.lower().replace(" ", "_") if col != "id" else col)
        for col in df.columns
    ]

    # merge with data on author and time
    data = data.merge(df, on=["id"], how="left")

    return data


def _calculate_mesh_balance(data: pd.DataFrame, mesh_terms_dict: dict) -> pd.DataFrame:
    """
    Calculate the balance of MeSH terms for each entry in the given DataFrame.
    Args:
        data (pd.DataFrame): Input DataFrame containing 'id' and 'mesh_terms' columns.
        mesh_terms_dict (dict): Dictionary mapping MeSH terms to their respective groups.
    Returns:
        pd.DataFrame: DataFrame with the balance of MeSH terms for each entry, merged
            with the original data.
    """
    df = data.copy()

    # transform mesh_terms by mapping to term_group
    df["mesh_terms"] = df["mesh_terms"].apply(
        lambda x: [y[0] for y in x] if x is not None else []
    )
    df["mesh_terms"] = df["mesh_terms"].apply(
        lambda x: [mesh_terms_dict.get(y, "") for y in x] if x is not None else []
    )

    # explode the mesh_terms column into multiple rows
    df = df.explode("mesh_terms")

    # group by author and mesh_terms and calculate the count
    df = df.groupby(["id", "mesh_terms"]).size().reset_index(name="count")

    # calculate the share of each mesh term
    df["share"] = df.apply(lambda row: row["count"], axis=1)

    # pivot the DataFrame to get one column for each mesh term
    df = df.pivot(index=["id"], columns="mesh_terms", values="share")

    # reset index
    df.reset_index(inplace=True)

    # change column names to be camel case, and prefix with "mesh_"
    df.columns = ["mesh_" + col if col != "id" else col for col in df.columns]

    # fill NaN values with 0
    df = df.fillna(0)

    # merge with data on author and time
    data = data.merge(df, on="id", how="left")

    return data


def _get_patent_citations(
    data: pd.DataFrame, patents_data: pd.DataFrame
) -> pd.DataFrame:
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
        patent_matches.groupby(["id"])
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
    data = data.merge(patent_matches_grouped, on="id", how="left")

    # fill count with 0
    data["patent_count"] = data["patent_count"].fillna(0)

    return data


def _get_pdb_activity(
    data: pd.DataFrame, pdb_submissions: pd.DataFrame
) -> pd.DataFrame:
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

    # extract from authorships every first item in the sublists
    pdb_submissions["authors"] = pdb_submissions["authorships"].apply(
        lambda x: [y[0] for y in x] if x is not None else []
    )

    # explode and drop id authors duplicates
    pdb_submissions_e = pdb_submissions.explode("authors")
    pdb_submissions_e.drop_duplicates(subset=["id", "authors"], inplace=True)

    # compute count, average resolution, and R_free per authors
    pdb_submissions_grouped_count = (
        pdb_submissions_e.groupby(["authors"])
        .agg(
            pdb_count=pd.NamedAgg(column="id", aggfunc="count"),
        )
        .reset_index()
    )

    # drop missing R_free and resolution values
    pdb_submissions_e["R_free"] = pd.to_numeric(
        pdb_submissions_e["R_free"], errors="coerce"
    )
    pdb_submissions_e["resolution"] = pd.to_numeric(
        pdb_submissions_e["resolution"], errors="coerce"
    )
    pdb_submissions_e = pdb_submissions_e[
        pdb_submissions_e["R_free"].notnull()
        & pdb_submissions_e["resolution"].notnull()
    ]

    # compute average resolution and R_free per authors
    pdb_submissions_grouped_chars = (
        pdb_submissions_e.groupby(["authors"])
        .agg(
            resolution=pd.NamedAgg(column="resolution", aggfunc="mean"),
            R_free=pd.NamedAgg(column="R_free", aggfunc="mean"),
        )
        .reset_index()
    )

    # merge the two groupeds
    pdb_submissions_grouped = pdb_submissions_grouped_count.merge(
        pdb_submissions_grouped_chars, on="authors", how="left"
    )

    # keep last author from data authorships
    data["last_author"] = data["authorships"].apply(
        lambda x: x[-1][0] if x is not None and len(x) > 0 else None
    )

    # keep a list of authors
    data_e = data.copy()
    data_e["authors"] = data_e["authorships"].apply(
        lambda x: [y[0] for y in x] if x is not None and len(x) > 0 else []
    )

    data_e = data_e.explode("authors")

    # merge with pdb_submissions_grouped
    data_e = data_e.merge(
        pdb_submissions_grouped, left_on="authors", right_on="authors", how="left"
    )

    # group by "id", sum counts, average R_free and resolution
    data_e_grouped = (
        data_e.groupby(["id"])
        .agg(
            pdb_count=pd.NamedAgg(column="pdb_count", aggfunc="sum"),
            resolution=pd.NamedAgg(column="resolution", aggfunc="mean"),
            R_free=pd.NamedAgg(column="R_free", aggfunc="mean"),
        )
        .reset_index()
    )

    # relabel to add "group_" to pdb_count, resolution, and R_free
    data_e_grouped.columns = [
        "group_" + col if col != "id" else col for col in data_e_grouped.columns
    ]

    data = data.merge(
        pdb_submissions[["id", "resolution", "R_free"]], on="id", how="left"
    )
    data = data.merge(data_e_grouped, on="id", how="left")

    return data
