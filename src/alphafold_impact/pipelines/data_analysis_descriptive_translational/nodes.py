"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
from joblib import Parallel, delayed
import altair as alt
from Bio import Entrez
from alphafold_impact.utils.altair import altair_to_png


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
    data = data[
        (data["level"] == "0") | (data["level"] == "1") | (data["level"] == "2")
    ]

    # if strength is not in the columns, add it as empty {}
    if "strength" not in data.columns:
        data["strength"] = "{}"

    data = data[
        [
            "id",
            "doi",
            "pmid",
            "parent_id",
            "publication_date",
            "parent_level",
            "level",
            "concepts",
            "mesh_terms",
            "topics",
            "strength",
        ]
    ]

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
    # preprocess data
    data = data.drop_duplicates(subset=["id", "parent_id", "level", "dataset"])
    data = data[["id", "doi", "pmid", "parent_id", "level", "dataset"]]

    # preprocess icite data
    icite_data["pmid"] = icite_data["pmid"].astype(str)
    icite_data = icite_data[["doi", "pmid", "cited_by_clin"]]

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

    def apply_func(row):
        return get_entrez_ptype_pmid(row)

    results = Parallel(n_jobs=8, verbose=10)(
        delayed(apply_func)(row) for row in combined_data["cited_by_clin"]
    )

    combined_data[["ca_publication_type", "ca_publication_date"]] = pd.DataFrame(
        results
    )

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


def get_patent_classifications(
    patent_data: pd.DataFrame,
    cpc_codes: pd.DataFrame,
):
    """
    Retrieves patent classifications from the given patent_data DataFrame and merges them
        with the cpc_codes DataFrame.

    Args:
        patent_data (pd.DataFrame): DataFrame containing patent data with columns
            'CPC Classifications' and 'IPCR Classifications'.
        cpc_codes (pd.DataFrame): DataFrame containing CPC codes and their corresponding
            titles.

    Returns:
        pd.DataFrame: DataFrame with patent classifications merged with cpc_codes, including
            additional columns 'parent_class_title' and 'class_title'.
    """

    # break CPC Classifications and IPCR Classifications based on ";;"
    patent_data["CPC Classifications"] = patent_data["CPC Classifications"].str.split(
        ";;"
    )
    patent_data["IPCR Classifications"] = patent_data["IPCR Classifications"].str.split(
        ";;"
    )

    # join them in a single list, removing duplicates
    patent_data["classifications"] = (
        patent_data["CPC Classifications"] + patent_data["IPCR Classifications"]
    )

    # drop the original columns
    patent_data.drop(
        columns=["CPC Classifications", "IPCR Classifications"], inplace=True
    )

    # explode the list
    patent_data = patent_data.explode("classifications")
    patent_data["classifications"] = patent_data["classifications"].astype(str)

    # drop parent_id, id, Title, classifications duplicates
    patent_data.drop_duplicates(
        subset=["parent_id", "id", "Title", "classifications"], inplace=True
    )

    # create column for parent classes (ie. before /)
    patent_data["parent_class"] = patent_data["classifications"].apply(lambda x: x[:4])

    # merge with cpc_codes
    cpc_codes["sort-key"] = cpc_codes["sort-key"].astype(str)
    patent_data = patent_data.merge(
        cpc_codes, left_on="parent_class", right_on="sort-key", how="left"
    )

    # rename columns
    patent_data.rename(columns={"text": "parent_class_title"}, inplace=True)
    patent_data.drop(columns="sort-key", inplace=True)

    # repeat for the full classification
    patent_data = patent_data.merge(
        cpc_codes, left_on="classifications", right_on="sort-key", how="left"
    )

    # rename columns
    patent_data.rename(columns={"text": "class_title"}, inplace=True)
    patent_data.drop(columns="sort-key", inplace=True)

    return patent_data


def create_tcc_sb_papers(
    data: pd.DataFrame,
):
    data["publication_date"] = pd.to_datetime(data["publication_date"])

    # Compute the number of days difference between 'publication_date' and 'ca_publication_date'
    data["tcc"] = (
        pd.to_datetime(data["ca_publication_date"]) - data["publication_date"]
    ).dt.days

    # change format of publication_date
    data["publication_date"] = data["publication_date"].dt.to_period("Q").astype(str)

    data = data[["publication_date", "source", "tcc"]]

    data = data[
        (data["publication_date"] >= "2021Q3") & (data["publication_date"] <= "2022Q4")
    ]
    chart = (
        (
            alt.Chart(data)
            .mark_boxplot()
            .encode(
                x=alt.X(
                    "source:O",
                    title=None,
                    axis=alt.Axis(labels=False, ticks=False),
                    scale=alt.Scale(padding=1),
                ),
                y=alt.Y(
                    "tcc:Q",
                    title="Median days to clinical citation",
                    scale=alt.Scale(zero=True),
                ),
                color=alt.Color(
                    "source:N",
                    title="Source",
                    # legend=alt.Legend(orient="top-right", offset=10),
                ),
                column=alt.Column(
                    "publication_date:N",
                    title="Source",
                ),
            )
        )
        .properties(width=100)
        .configure_facet(spacing=0)
        .configure_view(stroke=None)
    )

    chart = altair_to_png(chart)

    return data, chart


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


def create_publications_data(
    data: pd.DataFrame,
    source: str,
    mesh_terms: pd.DataFrame,
    patents_data,
    grants_data,
    pdb_submissions,
    icite_data,
):
    if "level" in data.columns:
        data = data[data["level"] != 3]
    try:
        data = _sort_drop(data)
    except:
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

    # grants NOTE: this should not be used for Other articles
    data = _get_awards(data, grants_data)

    # add pdb activity metrics
    data = _get_pdb_activity(data, pdb_submissions)

    # drop mesh_terms column
    data.drop(columns="mesh_terms", inplace=True)

    # icite data
    icite_outputs = create_cc_counts(
        data[["parent_id", "id", "level", "source", "pmid", "doi"]], icite_data
    )

    data = data.merge(icite_outputs[["id", "count"]], on="id", how="left")

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
    data_pmid = data_pmid.drop_duplicates(
        subset=["parent_id", "id", "level", "source"]
    )

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
    combined_data["count"] = combined_data["cited_by_clin"].apply(lambda x: len(x))

    # merge back with data, fill with 0 for missing count
    data = data.merge(combined_data[["id", "count"]], on="id", how="left")
    data["count"] = data["count"].fillna(0)

    return data


def merge_individual_data(data_af, data_ct, data_other):
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


def _get_field_distr(data):
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
        )
        .reset_index()
    )

    # merge with data
    data = data.merge(grant_matches_grouped, on=["id"], how="left")

    # fill count with 0
    data["grant_count"] = data["grant_count"].fillna(0)

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
