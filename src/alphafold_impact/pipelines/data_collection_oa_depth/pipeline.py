"""Pipeline for data collection from OpenAIRE API to a specific depth.

This pipeline fetches citation data from OpenAIRE API to a specific depth
and processes the data to create a network of citations.

To run this pipeline, use the following command:
    $ kedro run --pipeline=oa_data_collection_depth
"""

from kedro.pipeline import Pipeline, pipeline, node
from alphafold_impact import settings
from .nodes import (
    fetch_citation_all_depth,
    fetch_citation_to_specific_depth,
    process_data_by_level,
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

    mesh_processing = pipeline(
        [
            node(
                func=process_data_by_level,
                inputs={
                    "data": "raw_data",
                    "level": "level_placeholder",
                },
                outputs="intermediate_data",
            )
        ]
    )

    mesh_levels = [
        pipeline(
            mesh_processing,
            inputs={
                "raw_data": "oa.data_collection.depth.level.raw",
                "level_placeholder": f"params:oa.data_collection.depth.levels.{level}",
            },
            outputs={
                "intermediate_data": f"oa.data_collection.depth.mesh.{level}.intermediate",
            },
            tags=[
                f"oa.data_collection.depth.mesh.level.{str(level)}",
                "oa.data_collection.depth.mesh.levels",
            ],
            namespace=f"oa.data_collection.depth.mesh.level.{str(level)}",
        )
        for level in settings.DYNAMIC_PIPELINES_MAPPING["depth_levels"]
    ]

    return full_depth_pipeline + fixed_depth_pipeline + sum(mesh_levels)
