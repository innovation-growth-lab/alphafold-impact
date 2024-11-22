"""
This is a boilerplate pipeline 'data_analysis_chains'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from ..a_data_collection_oa.nodes import collect_papers  # pylint: disable=E0402

logger = logging.getLogger(__name__)


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

    logger.info(
        "Processed %d/%d papers with pmids from openalex",
        len(processed_papers_with_pmids),
        len(pmids),
    )

    # filter out from pdb_df the rows with pmid in processed_papers_with_pmids
    spdb_df = pdb_df[
        ~pdb_df["pmid"].isin(
            processed_papers_with_pmids["pmid"].astype(float).astype(str)
        )
    ]

    # turns out, missing dois still have a pdb entry doi,
    # not included in the API response!
    spdb_df["doi"] = spdb_df.apply(
        lambda row: (
            f"10.2210/pdb{row['rcsb_id'].lower()}/pdb"
            if row["doi"] in ["", "None", "nan"]
            else row["doi"]
        ),
        axis=1,
    )

    logger.info(
        "Attempting to process %d papers with dois from openalex",
        len(spdb_df),
    )

    # extract the list of doi from spdb_df
    dois = spdb_df["doi"].unique().tolist()

    processed_papers_with_dois = _get_papers(dois, api_config, label="doi")

    logger.info(
        "Processed %d/%d papers with dois from openalex",
        len(processed_papers_with_dois),
        len(dois),
    )

    # concatenate the two frames, drop duplicate ids
    oa_pdb_df = pd.concat(
        [processed_papers_with_pmids, processed_papers_with_dois]
    ).drop_duplicates(subset=["id"])

    pdb_with_pmid = pdb_df[~pdb_df["pmid"].isin(["", "None", "nan"])]
    pdb_with_pmid["pmid"] = pdb_with_pmid["pmid"].astype(float).astype(int).astype(str)

    logger.info("Merging openalex data with pdb data")
    oa_pdb_df_pmid = oa_pdb_df.merge(
        pdb_with_pmid[["rcsb_id", "pmid", "resolution", "R_free"]],
        how="right",
        on="pmid",
    )

    # remove rows with missing ids
    oa_pdb_df_pmid.dropna(subset="id", inplace=True)

    def _is_pdb_doi(doi):
        return doi is not None and doi.startswith("10.2210/pdb")

    # remove duplicate rows if one is the pdb entry
    oa_pdb_df_pmid["is_pdb_doi"] = oa_pdb_df_pmid["doi"].apply(_is_pdb_doi)
    oa_pdb_df_pmid = (
        oa_pdb_df_pmid.sort_values(by="is_pdb_doi")
        .drop_duplicates(subset="rcsb_id", keep="first")
        .drop(columns=["is_pdb_doi"])
    )

    logger.info("Merging openalex data with pdb data")
    oa_pdb_df_doi = oa_pdb_df.merge(
        pdb_df[["rcsb_id", "doi", "resolution", "R_free"]], how="right", on="doi"
    )

    logger.info("Concatenating openalex data with pdb data")
    combined_df = pd.concat([oa_pdb_df_pmid, oa_pdb_df_doi])  # .drop_duplicates("id")
    combined_df.dropna(subset="id", inplace=True)

    # remove again duplicate rows if one is the pdb entry
    combined_df["is_pdb_doi"] = combined_df["doi"].apply(_is_pdb_doi)
    combined_df = (
        combined_df.sort_values(by="is_pdb_doi")
        .drop_duplicates(subset="rcsb_id", keep="first")
        .drop(columns=["is_pdb_doi"])
    )

    logger.info("Finished processing %d papers", len(combined_df))
    # replace openalex.org in id
    combined_df["id"] = combined_df["id"].str.replace("https://openalex.org/", "")

    # drop ids
    combined_df.drop(columns=["ids"], inplace=True)

    return combined_df


def merge_uniprot_data(pdb_df: pd.DataFrame, uniprot_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merges Uniprot data with PDB data.

    Args:
        pdb_df (pd.DataFrame): The DataFrame containing PDB details.
        uniprot_df (pd.DataFrame): The DataFrame containing Uniprot details.

    Returns:
        pd.DataFrame: The combined DataFrame with PDB and Uniprot details.
    """

    logger.info("Transformations to PDB data")
    pdb_df["R_free"] = pd.to_numeric(pdb_df["R_free"], errors="coerce")
    pdb_df["resolution"] = pd.to_numeric(pdb_df["resolution"], errors="coerce")

    pdb_df = pdb_df.assign(
        publication_date=pd.to_datetime(pdb_df["publication_date"]),
        quarter=pd.PeriodIndex(pdb_df["publication_date"], freq="Q"),
    )

    pdb_oa_means = pdb_df.groupby("id").agg(
        resolution_mean=("resolution", "mean"),
        R_free_mean=("R_free", "mean"),
    )

    logger.info("merging to uniprot to get per-PDB primary status")
    merged_uniprot_df = uniprot_df.merge(
        pdb_df, how="left", left_on="pdb_id", right_on="rcsb_id"
    )

    # create column primary_submission for a given uniprot_id
    merged_uniprot_df["primary_submission"] = (
        merged_uniprot_df.groupby("uniprot_id")["publication_date"]
        .transform("min")
        .eq(merged_uniprot_df["publication_date"])
    )

    # create unique assignment
    primary_pdbs = (
        merged_uniprot_df[["rcsb_id", "primary_submission"]]
        .sort_values("primary_submission")
        .drop_duplicates(subset=["rcsb_id"], keep="last")
    )

    # transform boolean to int
    primary_pdbs["primary_submission"] = primary_pdbs["primary_submission"].astype(
        "Int64"
    )

    # merge back
    pdb_df = pdb_df.merge(primary_pdbs, how="left")

    # merge uniprot data with pdb data
    intermediate_df = uniprot_df.merge(
        pdb_df, how="left", left_on="pdb_id", right_on="rcsb_id"
    )

    logger.info("Creating disease counts in intermediate_df")
    intermediate_df["num_diseases"] = intermediate_df["disease_triples"].apply(
        lambda x: len(x) if x is not None else 0
    )

    logger.info("Creating inverse frequency of organism occurrence")
    organism_frequencies = (
        merged_uniprot_df.groupby(["organism_id", "quarter"])["organism_id"]
        .nunique()
        .groupby(level=0)
        .cumsum()
        .reset_index(name="cumulative_unique_organism_count")
    )
    organism_frequencies["organism_rarity"] = (
        1 / organism_frequencies["cumulative_unique_organism_count"]
    )

    # merge organism_rarity to intermediate_df
    intermediate_df = intermediate_df.merge(
        organism_frequencies[["organism_id", "quarter", "organism_rarity"]],
        how="left",
        on=["organism_id", "quarter"],
    )

    logger.info("Creating publication-level features")
    oa_structural_df = (
        intermediate_df.groupby("id")
        .agg(
            num_uniprot_structures=("uniprot_id", "nunique"),
            num_pdb_ids=("rcsb_id", "nunique"),
            num_primary_submissions=("primary_submission", "sum"),
            score_mean=("score", "mean"),
            complexity_sum=("complexity", "sum"),
            complexity_mean=("complexity", "mean"),
            organism_rarity_mean=("organism_rarity", "mean"),
            organism_rarity_max=("organism_rarity", "max"),
            num_diseases=("num_diseases", "sum"),
        )
        .reset_index()
    )

    # merge grouped_df with pdb_oa_means
    oa_structural_df = oa_structural_df.merge(
        pdb_oa_means, how="left", left_on="id", right_index=True
    )

    return oa_structural_df


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


def _get_papers(ids: list, api_config: dict, label: str = "doi"):
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
