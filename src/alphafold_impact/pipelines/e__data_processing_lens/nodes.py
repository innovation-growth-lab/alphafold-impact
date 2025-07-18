"""
This is a boilerplate pipeline 'data_analysis_exploratory'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict
import pandas as pd
from kedro.io import AbstractDataset

logger = logging.getLogger(__name__)


def create_and_explode_manual_data(
    data_loaders: Dict[str, AbstractDataset],
) -> pd.DataFrame:
    """
    Create and explode manual data by loading patents from data loaders, filtering and transforming the data.

    Args:
        data_loaders (Dict[str, AbstractDataset]): A dictionary of data loaders.

    Returns:
        pd.DataFrame: The concatenated and transformed dataframe containing all the patents.
    """

    # Initialize an empty list to store the dataframes
    dfs = []

    # Loop through each loader in patents
    for i, loader in enumerate(data_loaders.values()):

        logger.info("Loading patents from loader %s / %s", i + 1, len(data_loaders))

        # Call the loader function to get a dataframe
        df = loader()

        # Handle Spanish column names if present
        spanish_to_english = {
            "Jurisdicción": "Jurisdiction",
            "Tipo": "Kind",
            "Fecha de publicación": "Publication Date",
            "Título": "Title",
            "Recuento de patentes de citas": "Cites Patent Count",
            "Citado por recuento de patentes": "Cited by Patent Count",
            "NPL ID externos resueltos": "NPL Resolved External ID(s)",
            "Clasificaciones CPC": "CPC Classifications",
            "Clasificaciones IPCR": "IPCR Classifications",
            "Estado legal": "Legal Status",
        }

        # Rename Spanish columns if they exist
        df = df.rename(columns=spanish_to_english)

        # Keep the ones that have NPL Resolved External ID(s)
        df = df[df["NPL Resolved External ID(s)"].notnull()]

        # Keep only the necessary columns
        df = df[
            [
                "Jurisdiction",
                "Kind",
                "Publication Date",
                "Title",
                "Cites Patent Count",
                "Cited by Patent Count",
                "NPL Resolved External ID(s)",
                "CPC Classifications",
                "IPCR Classifications",
                "Legal Status",
            ]
        ]

        # Create a list from the NPL Resolved External ID(s) column, using as separator ;;
        df["NPL Resolved External ID(s)"] = df["NPL Resolved External ID(s)"].apply(
            lambda x: x.split(";;")
        )

        # Explode these lists
        df = df.explode("NPL Resolved External ID(s)")

        # Append the resulting dataframe to the list
        dfs.append(df)

    # Concatenate all the dataframes
    all_patents = pd.concat(dfs, ignore_index=True)
    return all_patents
