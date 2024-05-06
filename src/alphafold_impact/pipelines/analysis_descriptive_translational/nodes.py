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