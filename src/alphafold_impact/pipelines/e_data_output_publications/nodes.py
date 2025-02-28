"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from ..f_data_collection_foundational_labs.utils import (  # pylint: disable=E0402
    fetch_institution_data,
)

logger = logging.getLogger(__name__)


def update_alphafold_triad(data: pd.DataFrame) -> pd.DataFrame:
    """
    Uodates levels for AlphaFold triad other rows.

    Args:
        data (pd.DataFrame): The input DataFrame containing publication data.

    Returns:
        pd.DataFrame: The updated DataFrame with the new AlphaFold publication row
            appended and specified rows updated.
    """

    # change level to -1 for the other two papers
    data.loc[
        data["id"].isin(["W3211795435", "W3202105508", "W3177828909"]),
        ["level", "parent_level", "parent_id"],
    ] = [-1, -2, np.nan]

    return data


def create_publications_data(
    data: pd.DataFrame,
    source: str,
    mesh_terms: pd.DataFrame,
    patents_data: pd.DataFrame,
    pdb_submissions: pd.DataFrame,
    icite_data: pd.DataFrame,
    **kwargs,
) -> pd.DataFrame:
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

    if "seed_papers" in kwargs:
        logger.info("Assigning seed technologies")
        data = _assign_seed_technologies(data, kwargs["seed_papers"])

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
    if "seed_papers" in kwargs:
        data["source"] = source + "_" + data["ct_tech"]

    # create a new column 'parent_publication_date'
    id_date_dict = data.set_index("id")["publication_date"].to_dict()
    data["parent_publication_date"] = data["parent_id"].map(id_date_dict)

    # get short term citation counts
    # data = _normalise_citation_counts(data)

    # get fields data
    data = _get_field_distr(data)

    # get mesh terms
    data = _calculate_mesh_balance(data, mesh_terms_dict)

    # fill columns that begin with "mesh_" with 0
    data[[col for col in data.columns if col.startswith("mesh_")]] = data[
        [col for col in data.columns if col.startswith("mesh_")]
    ].fillna(0)

    # get patent citations
    data = get_patent_citations(data, patents_data)

    # add pdb activity metrics
    data = _get_pdb_activity(data, pdb_submissions)

    # drop mesh_terms column
    data.drop(columns="mesh_terms", inplace=True)

    # icite data
    icite_outputs = create_cc_counts(
        data[["parent_id", "id", "level", "source", "pmid", "doi"]], icite_data
    )

    data = data.merge(icite_outputs[["id", "ca_count"]], on="id", how="left")

    # id, parent id as str
    data["id"] = data["id"].astype(str)
    data["parent_id"] = data["parent_id"].astype(str)

    return data


def _sort_drop(data: pd.DataFrame) -> pd.DataFrame:
    """
    Sorts and filters a DataFrame based on specific criteria.
    Args:
        data (pd.DataFrame): The input DataFrame to be sorted and filtered.
    Returns:
        pd.DataFrame: The sorted and filtered DataFrame.
    """

    sort_order = {
        "strong": 0,
        "partial_strong": 1,
        "mixed": 2,
        "weak": 3,
        "partial_weak": 4,
        "unknown": 5,
        "no_data": 6,
    }
    data["sort_order"] = data.apply(
        lambda row: -1 if row["level"] == -1 else sort_order.get(row["chain_label"], 7),
        axis=1,
    )

    data = data.sort_values("sort_order").groupby(["id"]).first().reset_index()

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


def create_cc_counts(data, icite_data):
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
    combined_data["ca_count"] = combined_data["cited_by_clin"].apply(len)

    # merge back with data, fill with 0 for missing count
    data = data.merge(combined_data[["id", "ca_count"]], on="id", how="left")
    data["ca_count"] = data["ca_count"].fillna(0)

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

    data = pd.concat([data_af, data_ct], ignore_index=True)

    # drop id in data_other if it is in data_af or data_ct
    data_other = data_other[~data_other["id"].isin(data["id"])]

    # concatenate, ignore index
    data = pd.concat([data, data_other], ignore_index=True)

    # drop column '0' if it exists
    if "0" in data.columns:
        data.drop(columns="0", inplace=True)

    # change column level to str
    data["level"] = data["level"].astype(str)

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

    # add primary field
    data["primary_field"] = data["topics"].apply(
        lambda x: (
            x[0][5]
            if isinstance(x, np.ndarray)
            and len(x) > 0
            and isinstance(x[0], np.ndarray)
            and len(x[0]) > 0
            else None
        )
    )

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


def get_patent_citations(
    data: pd.DataFrame, patents_data: pd.DataFrame
) -> pd.DataFrame:
    """
    Get patent citations for the given data using multiple identifier types (DOI, PMID, PMCID).

    Args:
        data (pandas.DataFrame): The input data containing paper identifiers.
        patents_data (pandas.DataFrame): The patent data with NPL citations.

    Returns:
        pandas.DataFrame: The data with patent citations merged and aggregated.
    """

    # Create separate matches for each ID type
    matches = []

    # Match DOIs
    doi_matches = data.merge(
        patents_data, left_on="doi", right_on="NPL Resolved External ID(s)", how="inner"
    )
    matches.append(doi_matches)

    # Match PMIDs
    pmid_matches = data.merge(
        patents_data,
        left_on="pmid",
        right_on="NPL Resolved External ID(s)",
        how="inner",
    )
    matches.append(pmid_matches)

    # Combine all matches and remove duplicates
    patent_matches = pd.concat(matches, ignore_index=True)
    patent_matches = patent_matches.drop_duplicates(subset=["id", "Title"])

    # Group by paper ID and aggregate patent information
    patent_matches_grouped = (
        patent_matches.groupby(["id"])
        .agg(
            patent_count=pd.NamedAgg(column="Title", aggfunc="count"),
            CPCs=pd.NamedAgg(
                column="CPC Classifications",
                aggfunc=lambda x: ";;".join(map(str, x)),
            ),
            patent_citation=pd.NamedAgg(column="Cited by Patent Count", aggfunc="sum"),
        )
        .reset_index()
    )

    # Merge with original data and fill missing values
    data = data.merge(patent_matches_grouped, on="id", how="left")
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
    data = data.merge(
        pdb_submissions[
            [
                "id",
                "num_uniprot_structures",
                "num_pdb_ids",
                "num_primary_submissions",
                "score_mean",
                "complexity_sum",
                "complexity_mean",
                "organism_rarity_mean",
                "organism_rarity_max",
                "num_diseases",
                "resolution_mean",
                "R_free_mean",
                "mean_tmscore",
                "max_tmscore",
                "normalised_mean_tmscore",
                "normalised_max_tmscore",
            ]
        ],
        on="id",
        how="left",
        indicator=True,
    )
    data["pdb_submission"] = data["_merge"] == "both"
    data.drop(columns=["_merge"], inplace=True)

    return data


def select_regression_columns(data: pd.DataFrame, columns: list) -> pd.DataFrame:
    """
    Selects columns relevant for regression analysis from the given DataFrame.

    Args:
        data (pd.DataFrame): The input DataFrame containing various metrics and annotations.
        columns (list): A list of columns to keep in the DataFrame.

    Returns:
        pd.DataFrame: A DataFrame containing only the columns relevant for regression analysis.
    """

    # Drop columns that are not relevant for regression analysis
    regression_data = data.drop(columns=columns)

    return regression_data


def _assign_seed_technologies(
    data: pd.DataFrame, seed_papers: pd.DataFrame
) -> pd.DataFrame:
    """
    Assigns seed technologies to the data based on the seed papers.

    Args:
        data (pd.DataFrame): The input DataFrame containing the data.
        seed_papers (pd.DataFrame): The DataFrame containing the seed papers.

    Returns:
        pd.DataFrame: The DataFrame with the seed technologies assigned.
    """

    def create_subgroup(ct_data, ct_data_previous, level):
        return ct_data[
            (ct_data["level"] == level)
            & (ct_data["parent_id"].isin(ct_data_previous["id"]))
        ]

    # create seed papers
    ai_seeds = seed_papers[seed_papers["label"].str.contains("ai")]
    noai_seeds = seed_papers[~seed_papers["label"].str.contains("ai")]

    # initialise level 0 papers
    ai_l0 = data[(data["parent_id"].isin(ai_seeds["parent_id"])) & (data["level"] <= 0)]
    ai_l0["ct_tech"] = "ai"

    noai_l0 = data[
        (data["parent_id"].isin(noai_seeds["parent_id"])) & (data["level"] <= 0)
    ]
    noai_l0["ct_tech"] = "noai"

    levels = [ai_l0, noai_l0]
    for level in range(1, 3):
        ai_level = create_subgroup(data, levels[-2], level)
        ai_level["ct_tech"] = "ai"
        noai_level = create_subgroup(data, levels[-1], level)
        noai_level["ct_tech"] = "noai"
        levels.extend([ai_level, noai_level])

    # concatenate all levels and remove duplicates
    data = (
        pd.concat(levels)
        .sort_values(["ct_tech", "level"])
        .drop_duplicates(subset=["parent_id", "id", "level", "chain_label"])
        .reset_index(drop=True)
    )

    return data


def _normalise_citation_counts(data: pd.DataFrame) -> pd.DataFrame:
    """
    Normalises the citation counts in the given DataFrame.

    Args:
        data (pd.DataFrame): The input DataFrame containing the data.

    Returns:
        pd.DataFrame: A DataFrame with normalised citation counts.
    """

    data = data.drop_duplicates(subset=["id"])

    # explode the 'counts_by_year' column into multiple rows
    data_exploded = data.explode("counts_by_year")

    data_exploded.reset_index(drop=True, inplace=True)

    # separate the dictionary into two columns: 'year' and 'cited_by_count'
    data_exploded[["year", "cited_by_count"]] = pd.json_normalize(
        data_exploded["counts_by_year"]
    )

    # convert the 'publication_date' column to datetime and extract the year
    data_exploded["publication_date"] = pd.to_datetime(
        data_exploded["publication_date"]
    ).dt.year

    # calculate the 't' value
    data_exploded["t"] = data_exploded["year"] - data_exploded["publication_date"]

    # drop the rows with NaN values in the 't' column
    data_exploded.dropna(subset=["t"], inplace=True)

    data_exploded["t"] = data_exploded["t"].astype(int)

    # pivot the DataFrame to get the 't' columns
    t_columns = data_exploded.pivot_table(
        index="id", columns="t", values="cited_by_count", aggfunc="sum"
    )

    # remove negative columns
    t_columns = t_columns.loc[:, (t_columns.columns >= 0)]

    # Convert columns to Int64 type to accept NAs
    t_columns = t_columns.astype(pd.Int64Dtype())

    # renaming
    t_columns.columns = ["cit_" + str(int(i)) for i in t_columns.columns]

    # fill remaining NaN values in 't0' with 0
    t_columns["cit_0"] = t_columns["cit_0"].fillna(0)

    # merge the 't' DataFrame back to the original DataFrame
    data = data.merge(t_columns, left_on="id", right_index=True, how="left")

    return data


def get_institution_info(
    publications: pd.DataFrame,
) -> pd.DataFrame:
    """
    Fetches institution information for authors.

    Args:
        author_ids (pd.DataFrame): DataFrame containing author IDs and institution information.

    Returns:
        pd.DataFrame: DataFrame containing institution data.

    """
    # keep last author from data authorships
    publications["last_author_institution"] = publications["authorships"].apply(
        lambda x: x[-1][1] if x is not None and len(x) > 0 else None
    )

    institutions = publications["last_author_institution"].unique().tolist()

    logger.info("Fetching publications for %d authors", len(institutions))

    # slice to create parallel jobs that produce slices
    slices = [institutions[i : i + 150] for i in range(0, len(institutions), 150)]

    logger.info("Slicing authors into %d slices", len(slices))
    institution_data_list = Parallel(n_jobs=8, verbose=10)(
        delayed(fetch_institution_data)(institution)
        for slice_ in slices
        for institution in slice_
    )

    # filter out None values
    institution_data_list = [data for data in institution_data_list if data is not None]

    # convert the list of dictionaries to a DataFrame
    data = pd.DataFrame(institution_data_list)

    # merge back on author_ids
    publications = publications.merge(
        data,
        left_on="last_author_institution",
        right_on="institution",
        how="left",
        suffixes=("", "_institution"),
    )

    return publications



def define_high_pdb_authors(
    data: pd.DataFrame, pdb_submissions: pd.DataFrame
) -> pd.DataFrame:
    """
    Define high PDB authors based on the publication counts in the
    fourth quantile.

    Args:
        authors (list): The list of authors.
        pdb_submissions (pd.DataFrame): The DataFrame containing the PDB submissions.

    Returns:
        pd.DataFrame: The updated DataFrame with the 'high_pdb' column.
    """
    author_pdb_submissions = pdb_submissions.dropna(subset=["authorships"])

    author_pdb_submissions["author"] = author_pdb_submissions["authorships"].apply(
        lambda x: [y[0] for y in x] if x is not None else []
    )

    # explode and drop id authors duplicates
    author_pdb_submissions = author_pdb_submissions.explode("author")
    author_pdb_submissions.drop_duplicates(subset=["id", "author"], inplace=True)

    # do the same with the publications
    # keep last author from data authorships
    data["last_author"] = data["authorships"].apply(
        lambda x: x[-1][0] if x is not None and len(x) > 0 else None
    )

    data["author"] = data["authorships"].apply(
        lambda x: [y[0] for y in x] if x is not None else []
    )
    authors = data.explode("author")
    authors_list = authors["author"].unique()

    submissions = author_pdb_submissions.copy()
    submissions_pre2021 = submissions[
        submissions["publication_date"] < "2021-01-01"
    ]
    
    # keep only authors in the data
    submissions_pre2021 = submissions_pre2021[submissions_pre2021["author"].isin(authors_list)]

    # make sure each author is only counted once per id
    submissions_pre2021 = submissions_pre2021.drop_duplicates(subset=["id", "author"])

    # group by author and sum pdb_submission
    submissions_pre2021 = (
        submissions_pre2021.groupby("author")["num_pdb_ids"].sum().reset_index()
    )

    # create q4_pdb_pre2021 boolean column
    submissions_pre2021["q4_pdb_pre2021"] = submissions_pre2021["num_pdb_ids"].apply(
        lambda x: x >= submissions_pre2021["num_pdb_ids"].quantile(0.75)
    )

    # map author num_pdb ids
    authors["num_pdb_ids_pre2021"] = authors["author"].map(
        submissions_pre2021.set_index("author")["num_pdb_ids"]
    ).fillna(0).astype(int)
    authors["q4_pdb_pre2021"] = authors["author"].map(
        submissions_pre2021.set_index("author")["q4_pdb_pre2021"]
    ).fillna(False)

    # group by paper
    authors_grouped_in_paper = authors.groupby("id").agg(
        num_pdb_ids_pre2021=pd.NamedAgg(column="num_pdb_ids_pre2021", aggfunc="sum"),
        q4_pdb_pre2021_any=pd.NamedAgg(column="q4_pdb_pre2021", aggfunc="any")
    ).reset_index()

    # merge back with data
    data = data.merge(
        authors_grouped_in_paper[["id", "num_pdb_ids_pre2021", "q4_pdb_pre2021_any"]],
        on="id",
        how="left",
    )

    # fill missing values
    data["num_pdb_ids_pre2021"] = data["num_pdb_ids_pre2021"].fillna(0).astype(int)
    data["q4_pdb_pre2021_any"] = data["q4_pdb_pre2021_any"].fillna(False)

    return data