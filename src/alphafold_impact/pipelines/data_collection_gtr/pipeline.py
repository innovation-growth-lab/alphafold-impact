"""Pipeline for data collection.

This pipeline fetches data from the GtR API and preprocesses it into a
format that can be used by the rest of the project. To run this pipeline,
use the following command:

    $ kedro run --pipeline data_collection_gtr
"""
from kedro.pipeline import Pipeline, node, pipeline
from kedro.pipeline.modular_pipeline import pipeline as mpl

from .nodes import fetch_gtr_data, preprocess_data_to_df  # pylint: disable=E0401

template_pipeline = pipeline(
    [
        node(
            func=fetch_gtr_data,
            inputs=["params:param_requests", "params:endpoint_template"],
            outputs="raw_output_template",
        ),
        node(
            func=preprocess_data_to_df,
            inputs=["raw_output_template", "params:endpoint_template"],
            outputs="intermediate_output_template",
        ),
    ]
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=W0613
    """Pipeline for data collection.

    Returns:
        Pipeline: The data collection pipeline.
    """

    pipelines = list()
    for endpoint in ['projects', 'publications', 'organisations', 'funds']:
        # Create a pipeline for each endpoint
        pipeline_name = f'gtr-{endpoint}'
        pipe = mpl(
            pipe=template_pipeline,
            parameters={
                "endpoint_template": f"params:endpoints.{endpoint}",
            },
            outputs={
                "raw_output_template": f"gtr_raw_{endpoint}",
                "intermediate_output_template": f"gtr_intermediate_{endpoint}",
            },
            tags=pipeline_name,
        )
        # Add the pipeline to the list of pipelines
        pipelines.append(pipe)

    return sum(pipelines)


# [TODO] I would like to have a pipeline that can iterate over a number of possible endpoints for GtR data,
#       but I'm not sure how to do this with Kedro. I think I need to define a pipeline for each endpoint,
#       and then have a pipeline that calls each of these pipelines in turn. For now the best I've managed
#       is to have two separate pipelines, one for organisations and one for funds. I can use tags to
#       distinguish between the two. There's probably a better solution using modular pipelines:
#       https://docs.kedro.org/en/stable/nodes_and_pipelines/modular_pipelines.html
