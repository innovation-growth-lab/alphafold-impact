"""
This is a boilerplate pipeline 'data_analysis_chains'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from ..data_collection_oa.nodes import collect_papers  # pylint: disable=E0402

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


def _get_papers(ids: list, api_config: dict, label: str = "doi"):
    # slice of dois
    slice_keys = ["|".join(ids[i : i + 50]) for i in range(0, len(ids), 50)]

    # create parallel batches from slices
    batch_keys = [slice_keys[i : i + 100] for i in range(0, len(slice_keys), 100)]

    papers_with_ids = Parallel(n_jobs=8, backend="loky", verbose=10)(
        delayed(collect_papers)(
            oa_ids=batch_key,
            mailto=api_config["mailto"],
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


def collect_pdb_details(pdb_df: pd.DataFrame, api_config: dict) -> pd.DataFrame:
    """
    Collects PDB details from a DataFrame and API configuration.

    Args:
        pdb_df (pd.DataFrame): The DataFrame containing PDB details.
        api_config (dict): The API configuration.

    Returns:
        pd.DataFrame: The combined DataFrame with collected PDB details.
    """

    # get the set of unique pmid
    pmids = (
        pdb_df["pmid"]
        .replace(["", "None", "nan"], np.nan)
        .dropna()
        .astype(float)
        .astype(int)
        .astype(str)
        .unique()
        .tolist()
    )

    processed_papers_with_pmids = _get_papers(pmids, api_config, label="ids.pmid")

    # filter out from pdb_df the rows with pmid in processed_papers_with_pmids
    spdb_df = pdb_df[
        ~pdb_df["pmid"].isin(
            processed_papers_with_pmids["pmid"].astype(float).astype(str)
        )
    ]

    # extract the list of doi from spdb_df
    dois = (
        spdb_df["doi"].replace(["", "None", "nan"], np.nan).dropna().unique().tolist()
    )

    processed_papers_with_dois = _get_papers(dois, api_config, label="doi")

    # concatenate the two frames, drop duplicate ids
    oa_pdb_df = pd.concat(
        [processed_papers_with_pmids, processed_papers_with_dois]
    ).drop_duplicates(subset=["id"])

    pdb_with_pmid = pdb_df[~pdb_df["pmid"].isin(["", "None", "nan"])]
    pdb_with_pmid["pmid"] = pdb_with_pmid["pmid"].astype(float).astype(int).astype(str)

    # merge to oa_pdb_df the dataframe pdb_df, first try on pmid, then on doi
    oa_pdb_df_pmid = oa_pdb_df.merge(
        pdb_with_pmid[["rcsb_id", "pmid", "resolution", "R_free"]],
        how="left",
        on="pmid",
    )

    oa_pdb_df_doi = oa_pdb_df.merge(
        pdb_df[["rcsb_id", "doi", "resolution", "R_free"]], how="left", on="doi"
    )

    # Combine the results
    combined_df = pd.concat([oa_pdb_df_pmid, oa_pdb_df_doi]).drop_duplicates("id")

    # replace openalex.org in id
    combined_df["id"] = combined_df["id"].str.replace("https://openalex.org/", "")

    # drop ids
    combined_df.drop(columns=["ids"], inplace=True)

    return combined_df
