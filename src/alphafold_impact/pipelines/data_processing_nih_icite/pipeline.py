"""This module contains the pipeline for data processing for iCite data.

To run this pipeline, use the following command:

    $ kedro run --pipeline data_processing_nih_icite
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import filter_and_combine_icite


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=filter_and_combine_icite,
                inputs="nih.data_collection.icite.raw",
                outputs="nih.data_processing.icite.intermediate",
            ),
        ]
    )
