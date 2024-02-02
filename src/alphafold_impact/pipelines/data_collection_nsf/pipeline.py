"""
This is a boilerplate pipeline 'data_collection_nsf'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import fetch_nsf_data  # pylint: disable=E0401
from alphafold_impact import settings

def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=W0613
    """Pipeline for data collection.

    Returns:
        Pipeline: The data collection pipeline.
    """
    template_pipeline = pipeline(
        [
            node(
                func=fetch_nsf_data,
                inputs=["params:api_config", "params:variables", "params:digit"],
                outputs="raw",
                name="fetch_gtr_data",
            ),
        ]
    )

    pipelines = [
        pipeline(
            template_pipeline,
            parameters={"params:api_config": "nsf.data_collection.api"},
            namespace=f"nsf.data_collection.{label}",
            tags=[label, "nsf"],
        )
        for label in settings.DYNAMIC_PIPELINES_MAPPING["nsf"]
    ]
    return sum(pipelines)

# Note: "params:" get resolved with the namespace. Inputs need to specify them being parameters,
# otherwise they are assumed to be datasets. Their namespace is not prepended, ie. needs to define
# params:gtr.data_collection.{endpoint}.param_requests. Parameters argument assumes "params:" 
# prefix.