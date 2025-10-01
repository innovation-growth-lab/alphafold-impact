"""
Utility functions for data processing in the foundational labs pipeline.
"""

import logging
import pandas as pd
import numpy as np
from Bio import Entrez

Entrez.email = "david.ampudia@nesta.org.uk"

logger = logging.getLogger(__name__)


def get_usage(data: pd.DataFrame) -> pd.DataFrame:
    """Get the usage of frontier and AF methods"""
    data_subset = data[
        ["id", "authorships", "level", "publication_date", "source"]
    ].copy()

    # get the authors
    author_data = (
        data_subset.drop_duplicates(subset=["id"])
        .assign(
            authorships_parsed=lambda x: x["authorships"].apply(
                lambda y: (
                    [auth.split(",") for auth in y.split("|")]
                    if isinstance(y, str) and y != ""
                    else []
                )
            )
        )
        .drop(columns=["authorships"])
        .explode("authorships_parsed")
        .assign(
            author=lambda x: x["authorships_parsed"].apply(
                lambda y: y[0] if isinstance(y, list) and len(y) > 0 else None
            ),
            institution=lambda x: x["authorships_parsed"].apply(
                lambda y: y[1] if isinstance(y, list) and len(y) > 1 else None
            ),
        )
        .dropna(subset=["author"])
        .drop_duplicates(subset=["id", "author"])
        .reset_index(drop=True)
        .drop(columns=["authorships_parsed"])
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

    author_data["cumulative_counts"] = author_data.groupby(["author", "source"])[
        "counts"
    ].cumsum()

    logger.info("Pivoting table")
    author_data = author_data.pivot_table(
        index=["author", "quarter"],
        columns="source",
        values="cumulative_counts",
        fill_value=0,
    ).reset_index()

    # make source cols int
    for col in ["af", "ct_ai", "ct_pp", "ct_sb", "other"]:
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
    data_subset = data[
        ["id", "authorships", "publication_date", "chain_label", "source"]
    ].copy()

    # get the authors
    author_data = (
        data_subset.drop_duplicates(subset=["id"])
        .assign(
            authorships_parsed=lambda x: x["authorships"].apply(
                lambda y: (
                    [auth.split(",") for auth in y.split("|")]
                    if isinstance(y, str) and y != ""
                    else []
                )
            )
        )
        .drop(columns=["authorships"])
        .explode("authorships_parsed")
        .assign(
            author=lambda x: x["authorships_parsed"].apply(
                lambda y: y[0] if isinstance(y, list) and len(y) > 0 else None
            ),
            institution=lambda x: x["authorships_parsed"].apply(
                lambda y: y[1] if isinstance(y, list) and len(y) > 1 else None
            ),
        )
        .dropna(subset=["author"])
        .drop_duplicates(subset=["id", "author"])
        .reset_index(drop=True)
        .drop(columns=["authorships_parsed"])
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
    author_labels["ct_pp_with_intent"] = (
        author_labels["ct_pp_strong"]
        + author_labels["ct_pp_weak"]
        + author_labels["ct_pp_mixed"]
    )
    author_labels["ct_sb_with_intent"] = (
        author_labels["ct_sb_strong"]
        + author_labels["ct_sb_weak"]
        + author_labels["ct_sb_mixed"]
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
            resolution=pd.NamedAgg(column="resolution", aggfunc="mean"),
            R_free=pd.NamedAgg(column="R_free", aggfunc="mean"),
            max_tmscore=pd.NamedAgg(column="max_tmscore", aggfunc="max"),
            max_score=pd.NamedAgg(column="max_score", aggfunc="max"),
            max_fident=pd.NamedAgg(column="max_fident", aggfunc="max"),
            normalised_max_tmscore=pd.NamedAgg(
                column="normalised_max_tmscore", aggfunc="max"
            ),
            normalised_max_score=pd.NamedAgg(
                column="normalised_max_score", aggfunc="max"
            ),
            normalised_max_fident=pd.NamedAgg(
                column="normalised_max_fident", aggfunc="max"
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
    df = (
        df.groupby(["author", "quarter", "mesh_terms"]).size().reset_index(name="count")
    )

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
        lambda x: any("C524204448" in element for element in x)
    )

    # Group the filtered dataframe by pi_id and calculate the share of COVID concepts
    covid_share_2020 = (
        data_2020.groupby("author")["covid_2020"].sum()
        / data_2020.groupby("author")["covid_2020"].count()
    )

    # Convert the Series to a DataFrame
    covid_share_2020 = covid_share_2020.reset_index()

    # Rename the column
    covid_share_2020.columns = ["author", "covid_share_2020"]

    # Create a dictionary from the DataFrame
    covid_share_dict = dict(
        zip(covid_share_2020["author"], covid_share_2020["covid_share_2020"])
    )

    data["covid_share_2020"] = data["author"].map(covid_share_dict)

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
