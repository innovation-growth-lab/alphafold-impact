"""
This module contains utility functions for the data processing pipeline.
"""

import logging
from typing import Dict
import pandas as pd
from joblib import Parallel, delayed
from ..a_data_collection_oa.nodes import collect_papers  # pylint: disable=E0402

logger = logging.getLogger(__name__)


def _process_responses(responses):
    output = []

    for children_list in responses:

        json_data = [
            {
                k: v
                for k, v in item.items()
                if k
                in [
                    "id",
                    "ids",
                    "doi",
                    "publication_date",
                    "mesh_terms",
                    "cited_by_count",
                    "authorships",
                    "topics",
                    "concepts",
                    "fwci",
                ]
            }
            for item in children_list
        ]

        # transform to datafram
        df = pd.DataFrame(json_data)

        # if dataframe is empty, continue
        if df.empty:
            continue

        # extract pmid from ids
        df["pmid"] = df["ids"].apply(
            lambda x: (
                x.get("pmid").replace("https://pubmed.ncbi.nlm.nih.gov/", "")
                if x and x.get("pmid")
                else None
            )
        )

        # keep only a list of tuples with "descriptor_ui" and "descriptor_name" for each
        df["mesh_terms"] = df["mesh_terms"].apply(
            lambda x: (
                [(c["descriptor_ui"], c["descriptor_name"]) for c in x] if x else None
            )
        )

        # create boolean variable for neglected_disease if mesh_terms contains D058069
        df["neglected_disease"] = df["mesh_terms"].apply(
            lambda x: any([term[0] == "D058069" for term in x]) if x else False
        )

        # create a boolean variable for rare disease if mesh_terms contains D035583
        df["rare_disease"] = df["mesh_terms"].apply(
            lambda x: any([term[0] == "D035583" for term in x]) if x else False
        )

        # break atuhorship nested dictionary jsons, create triplets of authorship
        df["authorships"] = df["authorships"].apply(
            lambda x: (
                [
                    (
                        (
                            author["author"]["id"].replace("https://openalex.org/", ""),
                            inst["id"].replace("https://openalex.org/", ""),
                            author["author_position"],
                        )
                        if author["institutions"]
                        else [
                            author["author"]["id"].replace("https://openalex.org/", ""),
                            "",
                            author["author_position"],
                        ]
                    )
                    for author in x
                    for inst in author["institutions"] or [{}]
                ]
                if x
                else None
            )
        )

        # change doi to remove the url
        df["doi"] = df["doi"].str.replace("https://doi.org/", "")

        # create a list of topics, with id (replacing openalex.org),
        df["topics"] = df["topics"].apply(
            lambda x: (
                [
                    (
                        topic["id"].replace("https://openalex.org/", ""),
                        topic["display_name"],
                        topic["subfield"]["id"].replace("https://openalex.org/", ""),
                        topic["subfield"]["display_name"],
                        topic["field"]["id"].replace("https://openalex.org/", ""),
                        topic["field"]["display_name"],
                        topic["domain"]["id"].replace("https://openalex.org/", ""),
                        topic["domain"]["display_name"],
                    )
                    for topic in x
                ]
                if x
                else None
            )
        )

        # extract concepts, ie. for each element in the list of jsons
        df["concepts"] = df["concepts"].apply(
            lambda x: (
                [
                    (
                        concept["id"].replace("https://openalex.org/", ""),
                        concept["display_name"],
                    )
                    for concept in x
                ]
                if x
                else None
            )
        )

        # append to output
        output.append(df)

    df = pd.concat(output)

    return df


def get_papers(ids: list, api_config: dict, label: str = "doi") -> list:
    """
    Retrieve and process papers based on a list of IDs using a specified
      API configuration.

    Args:
        ids (list): A list of paper IDs to retrieve.
        api_config (dict): A dictionary containing API configuration
          parameters, including 'perpage'.
        label (str, optional): The filter criteria for the API request.
          Defaults to "doi".

    Returns:
        list: A list of processed papers with their corresponding IDs.
    """
    # slice of dois
    slice_keys = ["|".join(ids[i : i + 50]) for i in range(0, len(ids), 50)]

    # create parallel batches from slices
    batch_keys = [slice_keys[i : i + 100] for i in range(0, len(slice_keys), 100)]

    papers_with_ids = Parallel(n_jobs=8, backend="loky", verbose=10)(
        delayed(collect_papers)(
            oa_ids=batch_key,
            perpage=api_config["perpage"],
            filter_criteria=label,
            eager_loading=True,
            skip_preprocess=True,
        )
        for batch_key in batch_keys
    )

    # flatten the list of dicts into a single dict
    papers_with_ids_list = [v for d in papers_with_ids for _, v in d.items()]

    processed_papers_with_ids = _process_responses(papers_with_ids_list)

    return processed_papers_with_ids


def preprocess_pdb_dates(pdb_df: pd.DataFrame) -> Dict[str, pd.Timestamp]:
    """
    Preprocess PDB submission dates into a dictionary for quick lookup.
    """
    pdb_df["publication_date"] = pd.to_datetime(pdb_df["publication_date"])
    pdb_df["rcsb_id"] = pdb_df["rcsb_id"].str.lower()
    return dict(zip(pdb_df["rcsb_id"], pdb_df["publication_date"]))


def filter_and_compute_metrics(
    pairwise_similarities: pd.DataFrame, pdb_dates: Dict[str, pd.Timestamp]
) -> pd.DataFrame:
    """
    Filter similarity data for valid targets and compute metrics.

    Args:
        similarity_chunk (pd.DataFrame): The DataFrame containing similarity data.
        pdb_dates (Dict[str, pd.Timestamp]): A dictionary mapping PDB IDs to publication dates.

    Returns:
        pd.DataFrame: The DataFrame containing computed metrics.
    """
    pairwise_similarities["alntmscore"] = pairwise_similarities["alntmscore"].astype(
        float
    )
    pairwise_similarities["query_date"] = pairwise_similarities["query"].map(pdb_dates)
    pairwise_similarities["target_date"] = pairwise_similarities["target"].map(
        pdb_dates
    )

    # filter valid targets (target must be earlier or equal to query date, not self)
    valid_targets = pairwise_similarities[
        (pairwise_similarities["target_date"] < pairwise_similarities["query_date"])
        & (pairwise_similarities["query"] != pairwise_similarities["target"])
    ]

    # extract year from target_date
    valid_targets["query_year"] = valid_targets["query_date"].dt.year

    # compute yearly statistics (mean, std) for alntmscore
    yearly_stats = (
        valid_targets.groupby("query_year")["alntmscore"]
        .agg(["mean", "std"])
        .reset_index()
    )
    yearly_stats.rename(columns={"mean": "year_mean", "std": "year_std"}, inplace=True)

    # merge yearly stats back into valid_targets
    valid_targets = valid_targets.merge(yearly_stats, on="query_year", how="left")

    # compute normalised scores (Z-scores)
    valid_targets["normalised_score"] = (
        (valid_targets["alntmscore"] - valid_targets["year_mean"])
        / valid_targets["year_std"]
    ).fillna(
        0
    )  # Fill NaN for years with no variance

    # Compute metrics for each query
    def compute_metrics(group):
        top_25 = group.nlargest(25, "alntmscore")
        mean_tmscore = top_25["alntmscore"].mean()
        max_tmscore = group["alntmscore"].max()
        normalised_mean = top_25["normalised_score"].mean()
        normalised_max = group["normalised_score"].max()
        return pd.Series(
            {
                "mean_tmscore": mean_tmscore,
                "max_tmscore": max_tmscore,
                "normalised_mean_tmscore": normalised_mean,
                "normalised_max_tmscore": normalised_max,
            }
        )

    metrics = valid_targets.groupby("query").apply(compute_metrics).reset_index()

    # Ensure query IDs are uppercase
    metrics["query"] = metrics["query"].str.upper()

    return metrics
