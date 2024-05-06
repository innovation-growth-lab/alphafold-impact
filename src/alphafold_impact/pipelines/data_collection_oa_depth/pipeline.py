"""Pipeline for data collection from OpenAIRE API to a specific depth.

This pipeline fetches citation data from OpenAIRE API to a specific depth
and processes the data to create a network of citations.

To run this pipeline, use the following command:
    $ kedro run --pipeline data_collection_oa_depth
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    fetch_citation_all_depth,
    fetch_citation_to_specific_depth,
    preprocess_baseline_data,
    preprocess_restart_data,
    fetch_level_2_ct_papers
)


def create_pipeline(  # pylint: disable=C0116,W0613
    **kwargs,
) -> Pipeline:
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
                    "seed_paper": "params:oa.data_collection.depth.get.work_ids",
                    "api_config": "params:oa.data_collection.depth.api",
                    "filter_config": "params:oa.data_collection.depth.filter",
                    "max_depth": "params:oa.data_collection.depth.max_depth",
                },
                outputs="oa.data_collection.triad.depth.level.raw",
            )
        ],
        tags="oa.data_collection.depth.level",
    )

    fixed_depth_baseline_pipeline = pipeline(
        [
            node(
                func=preprocess_baseline_data,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.raw",
                    "alphafold_papers": "params:oa.data_collection.depth.get.work_ids",
                },
                outputs="oa.structural_biology.seed_papers",
            ),
            node(
                func=fetch_citation_to_specific_depth,
                inputs={
                    "seed_paper": "oa.structural_biology.seed_papers",
                    "api_config": "params:oa.data_collection.depth.api",
                    "filter_config": "params:oa.data_collection.depth.filter",
                    "max_depth": "params:sb_labs.data_collection.depth.max_depth",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.raw",
            ),
        ],
        tags="oa.data_collection.depth.structural_biology",
    )

    restart_depth_from_level_pipeline = pipeline(
        [
            node(
                func=preprocess_restart_data,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.raw",
                    "start_level": "params:start_level",
                },
                outputs=["seen_papers", "papers_to_process", "level"],
            ),
            node(
                func=fetch_citation_to_specific_depth,
                inputs={
                    "seed_paper": "papers_to_process",
                    "api_config": "params:oa.data_collection.depth.api",
                    "filter_config": "params:oa.data_collection.depth.filter",
                    "max_depth": "params:max_level",
                    "start_level": "level",
                    "papers_seen": "seen_papers",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.restarted.raw",
            ),
        ],
        tags="oa.data_processing.depth.structural_biology.restart",
    )

    fetch_level_2_ct_papers_pipeline = pipeline(
        [
            node(
                func=fetch_level_2_ct_papers,
                inputs={
                    "ct_data": "oa.data_processing.structural_biology.depth.reassigned.ct.intermediate",
                },
                outputs="ct_level_1_seeds",
            ),
            node(
                func=fetch_citation_to_specific_depth,
                inputs={
                    "seed_paper": "ct_level_1_seeds",
                    "start_level": "params:sb_labs.data_collection.depth.ct_start_level",
                    "api_config": "params:oa.data_collection.depth.api",
                    "filter_config": "params:oa.data_collection.depth.filter",
                    "max_depth": "params:sb_labs.data_collection.depth.max_depth",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.level.2.raw",
            ),
        ],
        tags="oa.data_collection.depth.level_2",
    )

    return full_depth_pipeline + fixed_depth_pipeline + fixed_depth_baseline_pipeline + restart_depth_from_level_pipeline + fetch_level_2_ct_papers_pipeline
