"""
This is a boilerplate pipeline 'data_processing_lens'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import create_and_explode_manual_data


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    return pipeline(
        [
            node(
                func=create_and_explode_manual_data,
                inputs={"data_loaders": "lens.data_processing.csv.raw"},
                outputs="lens.data_processing.primary",
                name="preprocess_lens_manual_data",
            ),
        ],
        tags=["data_processing_lens", "rerun"],
    )
