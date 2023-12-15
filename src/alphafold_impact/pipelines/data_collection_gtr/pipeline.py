"""Pipeline for data collection.

This pipeline fetches data from the GtR API and preprocesses it into a
format that can be used by the rest of the project. To run this pipeline,
use the following command:

    $ kedro run --pipeline data_collection_gtr
"""
from kedro.pipeline import Pipeline, node, pipeline

from .nodes import fetch_gtr_data, preprocess_data_to_df, preprocess_organisations  # pylint: disable=E0401


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=W0613
    """Pipeline for data collection.

    Returns:
        Pipeline: The data collection pipeline.
    """    
    return pipeline(
        [
            node(
                func=fetch_gtr_data,
                inputs="params:param_requests",
                outputs="gtr_raw_organisations",
                name="fetch_gtr_data_node",
            ),
            node(
                func=preprocess_data_to_df,
                inputs="gtr_raw_organisations",
                outputs="gtr_df_organisations",
                name="preprocess_data_to_df_node",
            ),
            node(
                func=preprocess_organisations,
                inputs="gtr_df_organisations",
                outputs="gtr_organisations",
                name="preprocess_organisations_node",
            ),
        ]
    )
