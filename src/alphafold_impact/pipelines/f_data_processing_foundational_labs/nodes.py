"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
from kedro.io import AbstractDataset
from .utils import (
    get_usage,
    get_intent,
    collect_covid_references,
    calculate_mesh_balance,
    process_institutional_data,
    process_pdb_data,
)
from ..g_data_collection_authors.nodes import ( # pylint: disable=E0402
    define_high_pdb_authors,
    calculate_field_share,
    get_patent_citations,
    create_cc_counts,
)

logger = logging.getLogger(__name__)


def get_lab_individual_outputs(
    data_loaders: AbstractDataset,
    publications_data: pd.DataFrame,
    pdb_submissions: pd.DataFrame,
    patents_data: pd.DataFrame,
    icite_data: pd.DataFrame,
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
    logger.info("Preparing data for end-of-loop merge")

    # process institutional data
    institutional_data = process_institutional_data(institutional_data)

    # create author submissions
    author_submissions = process_pdb_data(pdb_submissions)

    # drop publication_date from pdb_submissions
    pdb_submissions = pdb_submissions.drop(columns=["publication_date"])

    mesh_terms["term_group"] = mesh_terms["tree_number"].apply(
        lambda x: str(x)[:1] if x is not None else None
    )

    mesh_terms_dict = mesh_terms.set_index("DUI")["term_group"].to_dict()

    # create usage and intent data
    author_usage = get_usage(publications_data)
    author_intent = get_intent(publications_data)
    use_and_intent_cols_to_fill = ["af", "ct_ai", "ct_noai", "other"] + [
        col for col in author_intent.columns if col not in ["author", "quarter"]
    ]

    # extract author id (first item in each sublist in authorships)
    author_data = (
        publications_data[["authorships"]]
        .explode("authorships")
        .assign(
            author=lambda x: x["authorships"].apply(
                lambda y: y[0] if y is not None else None
            ),
            institution=lambda x: x["authorships"].apply(
                lambda y: y[1] if y is not None else None
            ),
        )
        .dropna(subset=["author"])
        .query("institution != ''")
        .drop_duplicates(subset=["author"])
        .reset_index(drop=True)
        .drop(columns=["authorships"])
    )

    output_data = pd.DataFrame()
    for i, loader in enumerate(data_loaders.values()):

        logger.info(
            "Processing data loader %d / %d",
            i + 1,
            len(data_loaders),
        )
        data = loader()

        data["publication_date"] = pd.to_datetime(data["publication_date"])
        data["quarter"] = data["publication_date"].dt.to_period("Q")

        data["pmcid"] = data["ids"].apply(
            lambda x: (
                x.get("pmcid", None).replace(
                    "https://www.ncbi.nlm.nih.gov/pmc/articles/", ""
                )
                if x.get("pmcid")
                else None
            )
        )

        # Help use functions designed for authors
        data.rename(columns={"pi_id": "author"}, inplace=True)

        # drop rows with publication_date older than 2017-06-01
        data = data[data["publication_date"] >= "2017-06-01"]

        # define whether author is a high pdb author
        data = define_high_pdb_authors(data, author_submissions)

        logger.info("Adding usage data")
        data = data.merge(author_usage, on=["author", "quarter"], how="left")

        logger.info("Adding intent data")
        data = data.merge(author_intent, on=["author", "quarter"], how="left")

        # sort and ffill
        data = data.sort_values(by=["author", "quarter"])
        for col in use_and_intent_cols_to_fill:
            data[col] = data.groupby("author")[col].ffill().fillna(0).astype(int)

        logger.info("Adding primary field")
        data["primary_field"] = data["topics"].apply(
            lambda x: x[0][5] if x is not None and len(x) > 0 else ""
        )

        logger.info("Calculating field share")
        data = calculate_field_share(data)

        logger.info("Calculating mesh balance")
        data = calculate_mesh_balance(data, mesh_terms_dict)

        logger.info("Collecting COVID references")
        data = collect_covid_references(data)

        for col in ["concepts", "mesh_terms", "grants", "topics", "ids", "authorships"]:
            try:
                data.drop(
                    columns=[col],
                    inplace=True,
                )
            except KeyError:
                pass

        data = data.merge(author_data, on=["author"], how="left")
        data = data.merge(institutional_data, on="institution", how="left")

        logger.info("Merging data with patent citations")
        data = get_patent_citations(data, patents_data)

        logger.info("Merging data with pdb submissions")
        data = data.merge(
            pdb_submissions, on="id", how="left", indicator=True
        )  # TO DOUBLECHECK
        data["pdb_submission"] = data["_merge"].apply(
            lambda x: True if x == "both" else False
        )
        data = data.drop(columns=["_merge"])

        logger.info("Merging data with clinical citations")
        icite_outputs = create_cc_counts(data[["id", "doi", "pmid"]], icite_data)
        data = data.merge(icite_outputs, on="id", how="left")

        # concatenate the data
        output_data = pd.concat([output_data, data], ignore_index=True)

    # create top quantile pre2021 pdb activity
    # get one row per author with their pre2021 pdb submissions
    author_pdb = output_data[
        ["author", "num_pdb_submissions_pre2021"]
    ].drop_duplicates()
    # calculate threshold based on author-level data
    pdb_threshold = author_pdb["num_pdb_submissions_pre2021"].quantile(0.75)
    # merge back to full dataset
    output_data["high_pdb_pre2021"] = (
        output_data["num_pdb_submissions_pre2021"] > pdb_threshold
    )

    return output_data
