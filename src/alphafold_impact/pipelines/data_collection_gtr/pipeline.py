"""Pipeline for data collection.

This pipeline fetches data from the GtR API and preprocesses it into a
format that can be used by the rest of the project. To run this pipeline,
use the following command:

    $ kedro run --pipeline data_collection_gtr

Alternatively, you can run this pipeline for a single endpoint:

    $ kedro run --pipeline data_collection_gtr --tags projects
"""
from kedro.pipeline import Pipeline, node, pipeline

from alphafold_impact import settings
from .nodes import fetch_gtr_data, preprocess_data_to_df  # pylint: disable=E0401


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=W0613
    """Pipeline for data collection.

    Returns:
        Pipeline: The data collection pipeline.
    """
    template_pipeline = pipeline(
        [
            node(
                func=fetch_gtr_data,
                inputs=["params_template", "endpoint_template"],
                outputs="raw",
                name="fetch_gtr_data",
            ),
            node(
                func=preprocess_data_to_df,
                inputs=["raw", "endpoint_template"],
                outputs="intermediate",
                name="preprocess_data_to_df",
            ),
        ]
    )

    pipelines = []
    for endpoint, label in settings.DYNAMIC_PIPELINES_MAPPING["gtr"]:
        pipelines.append(
            pipeline(
                template_pipeline,
                inputs={"params_template": f"params:gtr.{endpoint}.param_requests", "endpoint_template": f"params:gtr.{endpoint}.label"},
                namespace=f"gtr.{label}",
                tags=[endpoint, "gtr"],
            )
        )
    return sum(pipelines)
