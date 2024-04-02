"""Pipeline for data collection from OpenAIRE API to a specific depth.

This pipeline fetches citation data from OpenAIRE API to a specific depth
and processes the data to create a network of citations.

To run this pipeline, use the following command:
    $ kedro run --pipeline data_collection_oa_depth
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    fetch_citation_all_depth,
    fetch_citation_to_specific_depth
)


def create_pipeline(**kwargs) -> Pipeline:
    full_depth_pipeline = pipeline(
        [
            node(
                func=fetch_citation_all_depth,
                inputs={
                    "seed_paper": "params:get.work_id",
                    "api_config": "params:api",
                    "filter_config": "params:filter",
                },
                outputs=["edges", "works"],
                tags="network",
            )
        ],
        namespace="oa.data_collection.depth",
    )
    fixed_depth_pipeline = pipeline(
        [
            node(
                func=fetch_citation_to_specific_depth,
                inputs={
                    "seed_paper": "params:oa.data_collection.depth.get.work_id",
                    "api_config": "params:oa.data_collection.depth.api",
                    "filter_config": "params:oa.data_collection.depth.filter",
                    "max_depth": "params:oa.data_collection.depth.max_depth",
                },
                outputs="oa.data_collection.depth.level.raw",
            )
        ],
        tags="oa.data_collection.depth.level",
    )


    return (
        full_depth_pipeline
        + fixed_depth_pipeline
    )
