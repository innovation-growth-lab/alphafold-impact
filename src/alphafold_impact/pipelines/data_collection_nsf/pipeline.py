"""Pipeline for data collection.

This pipeline fetches data from the NSF API and preprocesses it into a
format that can be used by the rest of the project. To run this pipeline,
use the following command:

    $ kedro run --pipeline data_collection_nsf

Alternatively, you can run this pipeline for a single year:

    $ kedro run --pipeline data_collection_nsf --tags 2018
"""

from kedro.pipeline import Pipeline, pipeline, node
from alphafold_impact import settings
from .nodes import fetch_nsf_data  # pylint: disable=E0401


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=W0613
    """Pipeline for data collection.

    Returns:
        Pipeline: The data collection pipeline.
    """
    template_pipeline = pipeline(
        [
            node(
                func=fetch_nsf_data,
                inputs=["params:api_config", "params:variables", "raw"],
                outputs="intermediate",
                name="fetch_nsf_data",
            ),
        ]
    )

    pipelines = [
        pipeline(
            template_pipeline,
            parameters={
                "params:api_config": "nsf.data_collection.api",
                "params:variables": "nsf.data_collection.variables",
            },
            namespace=f"nsf.data_collection.awards.{label}",
            tags=[label, "nsf"],
        )
        for label in settings.DYNAMIC_PIPELINES_MAPPING["nsf"]
    ]
    return sum(pipelines)
