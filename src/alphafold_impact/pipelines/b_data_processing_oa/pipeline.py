"""This module contains the pipeline for data processing for OA data.

The pipeline contains nodes for processing data by level, combining data 
from different levels, and mesh tagging.

To run this pipeline, use the following command:

    $ kedro run --pipeline oa.data_processing.depth
"""

from kedro.pipeline import Pipeline, node, pipeline
from alphafold_impact import settings

from .nodes import (
    process_subfield_data,
    process_data_by_level,
    combine_levels_data,
    combine_levels_data_connect_parents,
    reassign_ct_levels,
    process_data_by_level_ptd,
    concat_pq_ptd,
)


def create_pipeline(  # pylint: disable=unused-argument&missing-function-docstring
    **kwargs,
) -> Pipeline:

    af_processing_pipeline = pipeline(
        [
            node(
                func=process_data_by_level,
                inputs={
                    "data": "oa.data_collection.triad.depth.level.raw",
                    "level": f"params:oa.data_collection.depth.levels.{level}",
                },
                outputs={
                    "intermediate_data": f"oa.data_processing.depth.no_mesh.{level}.intermediate",
                },
                tags=[
                    f"oa.data_processing.depth.no_mesh.level.{str(level)}",
                    "oa.data_processing.depth.no_mesh.levels",
                ],
                namespace=f"oa.data_processing.depth.no_mesh.level.{str(level)}",
            )
            for level in settings.DYNAMIC_PIPELINES_MAPPING["depth_levels"]
        ]
    )

    combine_levels_pipeline = pipeline(
        [
            node(
                func=combine_levels_data,
                inputs={
                    # "unique": "params:false_",
                    "level0": "oa.data_processing.depth.no_mesh.0.intermediate",
                    "level1": "oa.data_processing.depth.no_mesh.1.intermediate",
                    "level2": "oa.data_processing.depth.no_mesh.2.intermediate",
                },
                outputs="oa.data_processing.depth.intermediate",
            )
        ],
        tags="oa.data_processing.depth.combine_levels",
    )

    structural_biology_processing_pieline = pipeline(
        [
            node(
                func=process_subfield_data,
                inputs=[
                    "oa.data_collection.subfield.structural_biology.raw"
                ],
                outputs=[
                    "oa.data_processing.subfield.structural_biology.primary"
                ],
                tags=[
                    "oa.data_processing.subfield.structural_biology",
                ],
            )
        ]
    )

    baseline_level0_pipeline = pipeline(
        [
            node(
                func=process_data_by_level,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.raw",
                    "level": "params:oa.data_collection.depth.levels.0",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.0.intermediate",
            )
        ],
        tags=[
            "oa.data_processing.structural_biology.depth.level.0",
            "oa.data_processing.structural_biology.depth.levels",
        ],
    )

    baseline_level1_pipeline = pipeline(
        [
            node(
                func=process_data_by_level_ptd,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.raw",
                    "level": "params:oa.data_collection.depth.levels.1",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.1.ptd.intermediate",
            ),
            node(
                func=concat_pq_ptd,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.1.ptd.intermediate",  # pylint: disable=line-too-long
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.1.intermediate",
                tags="concat_pq_ptd",
            ),
        ],
        tags=[
            "oa.data_processing.structural_biology.depth.level.1",
            "oa.data_processing.structural_biology.depth.levels",
        ],
    )

    reassign_baseline_levels_pipeline = pipeline(
        [
            node(
                func=combine_levels_data,
                inputs={
                    "level_seed": "oa.data_processing.subfield.structural_biology.primary",
                    "level0": "oa.data_collection.subfield.structural_biology.depth.0.intermediate",
                    "level1": "oa.data_collection.subfield.structural_biology.depth.1.intermediate",
                },
                outputs="oa.data_processing.structural_biology.depth.intermediate",
            ),
            node(
                func=reassign_ct_levels,
                inputs={
                    "data": "oa.data_processing.structural_biology.depth.intermediate",
                    "ct_data": "chains.seed_technologies.intermediate",
                },
                outputs=[
                    "oa.data_processing.structural_biology.depth.reassigned.ct.intermediate",
                    "oa.data_processing.structural_biology.depth.other.intermediate",
                ],
            ),
        ],
        tags="oa.data_processing.depth.reassign_ct",
    )

    post_level2_dw_pipeline = pipeline(
        [
            node(
                func=process_data_by_level_ptd,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.raw",
                    "level": "params:oa.data_collection.depth.levels.2",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.2.ptd.intermediate",
            ),
            node(
                func=concat_pq_ptd,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.2.ptd.intermediate", # pylint: disable=line-too-long
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.2.intermediate",
                tags=["concat_pq_ptd", "concat_and_combine_ct"],
            ),
            node(
                func=combine_levels_data_connect_parents,
                inputs={
                    "level1": "oa.data_processing.structural_biology.depth.reassigned.ct.intermediate", # pylint: disable=line-too-long
                    "level2": "oa.data_collection.subfield.structural_biology.depth.2.intermediate",
                },
                outputs="oa.data_processing.structural_biology.depth.ct.intermediate",
                tags=["combine_ct", "concat_and_combine_ct"],
            ),
        ],
        tags="oa.data_processing.depth.post_level2_dw",
    )

    return (
        af_processing_pipeline
        + combine_levels_pipeline
        + structural_biology_processing_pieline
        + baseline_level0_pipeline
        + baseline_level1_pipeline
        + reassign_baseline_levels_pipeline
        + post_level2_dw_pipeline
    )
