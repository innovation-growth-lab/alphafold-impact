"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""
import logging
import pandas as pd

logger = logging.getLogger(__name__)

def load_input_data(
    data: pd.DataFrame,
    source: str,
):
    data["level"] = data["level"].astype(str)

    logger.info("Filter for level 0")
    data = data[data["level"] == "0"]

    logger.info("Create columns with source")
    data["source"] = source

    return data

def merge_inputs(**kwargs):
    return pd.concat(
        [kwargs["alphafold_data"], kwargs["ct_data"], kwargs["other_data"]]
    )

def get_cc_counts(data: pd.DataFrame, icite_data: pd.DataFrame):
    # change level to -1 for parent_id "W3177828909", "W3211795435", "W3202105508"
    data.loc[data["id"].isin(["W3177828909", "W3211795435", "W3202105508"]), "level"] = -1

    # preprocess icite data
    icite_data["pmid"] = icite_data["pmid"].astype(str)

    logger.info("Merging chains with iCite data")

    # Merge on 'pmid'
    data_pmid = data.merge(
        icite_data[['pmid', 'cited_by_clin']],
        how='left',
        on='pmid'
    )

    # Remove duplicates
    data_pmid = data_pmid.drop_duplicates(subset=["parent_id", "id", "level", "source"])

    # Merge on 'doi'
    data_doi = data.merge(
        icite_data[['doi', 'cited_by_clin']],
        how='left',
        on='doi'
    )

    # Remove duplicates
    data_doi = data_doi.drop_duplicates(subset=["parent_id", "id", "level", "source"])

    # Combine the results
    combined_data = pd.concat([data_pmid, data_doi]).drop_duplicates()

    # how many CC
        # strong and weak
    # how many strong articles
    # patenst

    return combined_data