"""
This is a boilerplate pipeline 'data_analysis_chains'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
import numpy as np
from scipy.spatial.distance import euclidean
from ..d_data_processing_chains.utils import (  # pylint: disable=relative-beyond-top-level
    sort_and_drop,
)


logger = logging.getLogger(__name__)


def _assign_label(row):
    """
    Assigns a label to a row based on the values of the 'parent_pp_topic',
    'parent_protein_concept', and 'parent_ai_concept' columns.

    Args:
        row (pandas.Series): The row containing the data.

    Returns:
        str: The assigned label based on the conditions.

    """
    if row["parent_pp_topic"] or row["parent_protein_concept"]:
        if not row["parent_ai_concept"]:
            return "pp"
        else:
            return "pp+ai"
    elif row["parent_ai_concept"]:
        return "sb+ai"
    else:
        return "sb"


def process_baseline_data(
    alphafold_data: pd.DataFrame,
    intent_data: pd.DataFrame,
    baseline_data: pd.DataFrame,
    seed_baseline_data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Process baseline data by merging it with seed data, grouping and aggregating the data,
    and returning the baseline candidates.

    Args:
        intent_data (pd.DataFrame): DataFrame containing intent data.
        baseline_data (pd.DataFrame): DataFrame containing baseline data.
        seed_baseline_data (pd.DataFrame): DataFrame containing seed baseline data.

    Returns:
        pd.DataFrame: DataFrame containing the baseline candidates.
    """
    baseline_data["parent_level"] = baseline_data["level"] - 1

    logger.info("Merging baseline data with seed baseline data")
    df = baseline_data.merge(
        seed_baseline_data[
            ["id", "doi", "pmid", "concepts", "topics", "publication_date"]
        ],
        left_on="parent_id",
        right_on="id",
        how="left",
        suffixes=("", "_parent"),
    )

    df.rename(
        columns={
            "doi_parent": "parent_doi",
            "pmid_parent": "parent_pmid",
            "concepts_parent": "parent_concepts",
            "topics_parent": "parent_topics",
            "publication_date_parent": "parent_publication_date",
        },
        inplace=True,
    )
    df.drop(columns=["id_parent"], inplace=True)

    logger.info("Creating dictionaries for doi matching")
    processed_grouped_data_doi = (
        intent_data.groupby(["parent_doi", "doi"])[["influential", "intent", "context"]]
        .apply(lambda x: x.to_dict(orient="records"))
        .reset_index()
        .rename(columns={0: "strength"})
    )

    logger.info("Creating dictionaries for pmid matching")
    processed_grouped_data_pmid = (
        intent_data.groupby(["parent_doi", "pmid"])[
            ["influential", "intent", "context"]
        ]
        .apply(lambda x: x.to_dict(orient="records"))
        .reset_index()
        .rename(columns={0: "strength"})
    )

    logger.info("Merging citation links with baseline data")
    processed_data = pd.merge(
        df, processed_grouped_data_doi, on=["parent_doi", "doi"], how="left"
    )

    processed_data = pd.merge(
        processed_data,
        processed_grouped_data_pmid,
        left_on=["parent_doi", "parent_pmid"],
        right_on=["parent_doi", "pmid"],
        how="left",
    )

    logger.info("Combining citation links")
    processed_data["strength"] = processed_data["strength_x"].combine_first(
        processed_data["strength_y"]
    )

    # rename the pmid columns
    processed_data.rename(columns={"pmid_x": "pmid"}, inplace=True)
    processed_data.drop(columns=["strength_x", "strength_y", "pmid_y"], inplace=True)

    # drop if id is W3177828909, W3211795435, W3202105508
    processed_data = processed_data[
        ~processed_data["id"].isin(["W3177828909", "W3211795435", "W3202105508"])
    ]

    logger.info("Exploding citation links")
    processed_data = processed_data.explode("strength")

    logger.info("Remove citation links with no strength")
    processed_data.dropna(subset=["strength"], inplace=True)

    for col in ["intent", "context"]:
        logger.info("Extracting %s from citation links", col)
        processed_data[col] = processed_data["strength"].apply(
            lambda x: x[col] if x else None  # pylint: disable=cell-var-from-loop
        )

    logger.info("Creating processing variables")
    processed_data["intent"] = processed_data["intent"].apply(
        lambda x: (
            "result"
            if x == "result"
            else (
                "methodology"
                if x == "methodology"
                else (
                    "background"
                    if x == "background"
                    else ("unknown" if x == "unknown" else "no_data")
                )
            )
        )
    )

    # define & implement the custom sorting, filtering order
    processed_data = sort_and_drop(processed_data, unique=False)

    # drop empty topics & concepts rows
    processed_data = processed_data.dropna(
        subset=["topics", "concepts", "parent_topics", "parent_concepts"]
    )

    # create a boolean column for whether topics includes "T10044", or "T12254"
    processed_data["parent_pp_topic"] = processed_data["parent_topics"].apply(
        lambda x: (
            any("T10044" in sublist[0] or "T12254" in sublist[0] for sublist in x)
        )
    )

    # create a boolean if it includes protein concepts, C18051474 and C47701112
    processed_data["parent_protein_concept"] = processed_data["parent_concepts"].apply(
        lambda x: (
            any("C18051474" in sublist[0] or "C47701112" in sublist[0] for sublist in x)
        )
    )

    # create a boolean if AI concept
    processed_data["parent_ai_concept"] = processed_data["parent_concepts"].apply(
        lambda x: any("C154945302" in sublist[0] for sublist in x)
    )

    # select only if parent_publication_date >= 2018-01-01 and <= 2022-06-15
    processed_data = processed_data[
        (processed_data["parent_publication_date"] >= "2018-01-01") &
        (processed_data["parent_publication_date"] <= "2022-06-15")
    ]

    logger.info("Creating aggregated baseline data")
    baseline_agg = (
        processed_data.groupby("parent_id")["intent"]
        .value_counts(normalize=True)
        .unstack()
        .fillna(0)
        .reset_index()
    )

    baseline_agg["num_citations"] = (
        processed_data.groupby("parent_id")["intent"].count().reset_index()["intent"]
    )

    # add parent_pp_topic, parent_protein_concept, parent_ai_concept
    baseline_agg = baseline_agg.merge(
        processed_data.groupby("parent_id")[
            ["parent_pp_topic", "parent_protein_concept", "parent_ai_concept"]
        ]
        .first()
        .reset_index(),
        on="parent_id",
    )

    # get candidates: more than 50 num_citations
    baseline_candidates = baseline_agg[baseline_agg["num_citations"] > 50]

    # manually remove past matches deemed not suitable W2973523639, W2984761660, W2999044305, W3104537585, W3110645309, W3199799076
    baseline_candidates = baseline_candidates[
        ~baseline_candidates["parent_id"].isin(
            [
                "W2973523639",
                "W2984761660",
                "W2999044305",
                "W3104537585",
                "W3110645309",
                "W3199799076",
                "W2944959599",
            ]
        )
    ]

    # remove candidates if parent_id is in "id" in alphafold_data
    baseline_candidates = baseline_candidates[
        ~baseline_candidates["parent_id"].isin(alphafold_data["id"])
    ]

    return baseline_candidates


def process_af_data(alphafold_data: pd.DataFrame) -> pd.DataFrame:
    """
    Process AlphaFold data by extracting information from citation links and
        calculating proportions.

    Args:
        alphafold_data (pd.DataFrame): The input DataFrame containing AlphaFold
            data.

    Returns:
        pd.DataFrame: The processed DataFrame with extracted information and
            calculated proportions.
    """

    # compare to alphafold baselines
    alphafold_df = alphafold_data[alphafold_data["level"] == 0]
    processed_alphafold = alphafold_df.explode("strength")

    # logger.info("Remove citation links with no strength")
    processed_alphafold.dropna(subset=["strength"], inplace=True)

    for col in ["intent", "context"]:
        logger.info("Extracting %s from citation links", col)
        processed_alphafold[col] = processed_alphafold["strength"].apply(
            lambda x: x[col] if x else None  # pylint: disable=cell-var-from-loop
        )

    # create strong variable
    processed_alphafold["intent"] = processed_alphafold["intent"].apply(
        lambda x: (
            "result"
            if x == "result"
            else (
                "methodology"
                if x == "methodology"
                else (
                    "background"
                    if x == "background"
                    else ("unknown" if x == "unknown" else "no_data")
                )
            )
        )
    )

    # get for each parent_id in processed_alphafold, the proportion of intents
    alphafold_target = (
        processed_alphafold.groupby("parent_id")["intent"]
        .value_counts(normalize=True)
        .unstack()
        .fillna(0)
        .reset_index()
    )

    # also add a column for number of citations for each parent_id
    alphafold_target["num_citations"] = (
        processed_alphafold.groupby("parent_id")["intent"]
        .count()
        .reset_index()["intent"]
    )

    return alphafold_target


def assign_focal_label(
    baseline_candidates: pd.DataFrame, alphafold_target: pd.DataFrame
) -> pd.DataFrame:
    """
    Assigns labels to baseline candidates based on the similarity of their results
    with the results from the AlphaFold model.

    Args:
        baseline_candidates (pd.DataFrame): DataFrame containing baseline candidates.
        alphafold_target (pd.DataFrame): DataFrame containing AlphaFold target results.

    Returns:
        pd.DataFrame: DataFrame with assigned labels.

    """
    distances = []
    for _, row1 in baseline_candidates.iterrows():
        row_distance = []
        for _, row2 in alphafold_target.iterrows():
            baseline_r = np.array([row1["result"], row1["methodology"]])
            alphafold_r = np.array([row2["result"], row2["methodology"]])
            distance = euclidean(baseline_r, alphafold_r)
            row_distance.append(distance)
        distances.append(row_distance)

    # assign the average of the inner lists back to baseline_candidates
    baseline_candidates["distance"] = np.mean(distances, axis=1)

    # drop any row with a distance larger than the quantile of the distances
    threshold = baseline_candidates["distance"].quantile(0.35)
    baseline_candidates_thresholded = baseline_candidates[
        baseline_candidates["distance"] < threshold
    ]

    baseline_candidates_thresholded["label"] = baseline_candidates_thresholded.apply(
        _assign_label, axis=1
    )

    return baseline_candidates_thresholded
