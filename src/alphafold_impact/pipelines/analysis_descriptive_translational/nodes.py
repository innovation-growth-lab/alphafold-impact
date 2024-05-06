"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""
import logging
import pandas as pd

logger = logging.getLogger(__name__)

def load_input_data(
    alphafold_data: pd.DataFrame,
    ct_data: pd.DataFrame,
    other_data: pd.DataFrame,
):
    
    logger.info("Filtering data by level 0")
    af_sb_data = alphafold_data[alphafold_data["level"]==0]
    ct_sb_data = ct_data[ct_data["level"]=="0"]
    other_sb_data = other_data[other_data["level"]=="0"]

    logger.info("Create columns with source")
    af_sb_data["source"] = "alphafold"
    ct_sb_data["source"] = "ct"
    other_sb_data["source"] = "other"

    logger.info("Concatenating data")
    level0_data = pd.concat([af_sb_data, ct_sb_data, other_sb_data])

    return level0_data