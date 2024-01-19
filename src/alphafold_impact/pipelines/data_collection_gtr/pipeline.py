"""Pipeline for data collection.

This pipeline fetches data from the GtR API and preprocesses it into a
format that can be used by the rest of the project. To run this pipeline,
use the following command:

    $ kedro run --pipeline data_collection_gtr

Alternatively, you can run this pipeline for a single endpoint:

    $ kedro run --pipeline data_collection_gtr --tags projects

In regards to the use of namespaces, note that these are appended as
prefixes to the outputs of the nodes in the pipeline. 

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
                inputs=["params:param_requests", "params:label", "params:test"],
                outputs="raw",
                name="fetch_gtr_data",
            ),
            node(
                func=preprocess_data_to_df,
                inputs=["raw", "params:label"],
                outputs="intermediate",
                name="preprocess_data_to_df",
            ),
        ]
    )

    pipelines = []
    for label in settings.DYNAMIC_PIPELINES_MAPPING["gtr"]:
        pipelines.append(
            pipeline(
                template_pipeline,
                parameters={"params:test": "test"},
                namespace=f"gtr.data_collection.{label}",
                tags=[label.rsplit('/', maxsplit=1)[-1], "gtr"],
            )
        )
    return sum(pipelines)

# Note: "params:" get resolved with the namespace. Inputs need to specify them being parameters,
# otherwise they are assumed to be datasets. Their namespace is not prepended, ie. needs to define
# params:gtr.data_collection.{endpoint}.param_requests. Parameters argument assumes "params:" 
# prefix.