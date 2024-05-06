"""This module contains the pipeline for data processing for OA data.

The pipeline contains nodes for processing data by level, combining data 
from different levels, and mesh tagging.

To run this pipeline, use the following command:

    $ kedro run --pipeline oa.data_processing.depth
"""

from kedro.pipeline import Pipeline, node, pipeline
from alphafold_impact import settings
from ...utils.nih_mesh_tagging import (  # pylint: disable=E0402
    generate_text_for_mesh_tagging,
    skr_web_python_api_generic_batch,
    mesh_results_to_df,
)
from .nodes import (
    process_subfield_data,
    process_data_by_level,
    combine_levels_data,
    reassign_ct_levels,
    process_data_by_level_ptd,
    concat_pq_ptd
)

def create_pipeline(  # pylint: disable=unused-argument&missing-function-docstring
    **kwargs,
) -> Pipeline:

    mesh_processing = pipeline(
        [
            node(
                func=process_data_by_level,
                inputs={
                    "data": "raw_data",
                    "level": "level_placeholder",
                    "extra_mesh": "extra_mesh_placeholder",
                },
                outputs="intermediate_data",
            )
        ]
    )

    additioal_mesh_levels = [
        pipeline(
            mesh_processing,
            inputs={
                "raw_data": "oa.data_collection.depth.level.raw",
                "level_placeholder": f"params:oa.data_collection.depth.levels.{level}",
                "extra_mesh_placeholder": "params:true_",
            },
            outputs={
                "intermediate_data": f"oa.data_processing.depth.mesh.{level}.intermediate",
            },
            tags=[
                f"oa.data_processing.depth.mesh.level.{str(level)}",
                "oa.data_processing.depth.mesh.levels",
            ],
            namespace=f"oa.data_processing.depth.mesh.level.{str(level)}",
        )
        for level in settings.DYNAMIC_PIPELINES_MAPPING["depth_levels"]
    ]

    no_additional_mesh_levels = [
        pipeline(
            mesh_processing,
            inputs={
                "raw_data": "oa.data_collection.triad.depth.level.raw",
                "level_placeholder": f"params:oa.data_collection.depth.levels.{level}",
                "extra_mesh_placeholder": "params:false_",
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

    combine_levels_pipeline = pipeline(
        [
            node(
                func=combine_levels_data,
                inputs={
                    # "unique": "params:false_",
                    "level0": "oa.data_processing.depth.no_mesh.0.intermediate",
                    "level1": "oa.data_processing.depth.no_mesh.1.intermediate",
                    "level2": "oa.data_processing.depth.no_mesh.2.intermediate",
                    "level3": "oa.data_processing.depth.no_mesh.3.intermediate",
                },
                outputs="oa.data_processing.depth.intermediate",
            )
        ],
        tags="oa.data_processing.depth.combine_levels",
    )

    mesh_tagging_pipeline = pipeline(
        [
            node(
                func=generate_text_for_mesh_tagging,
                inputs={
                    "df": "oa.data_processing.depth.intermediate",
                    "id_col": "params:oa.data_processing.depth.mesh_tagging.id_col",
                    "text_col": "params:oa.data_processing.depth.mesh_tagging.text_col",
                },
                outputs="oa.data_processing.depth.mesh_tagging.input",
            ),
            node(
                func=skr_web_python_api_generic_batch,
                inputs={
                    "text_to_tag": "oa.data_processing.depth.mesh_tagging.input",
                    "email": "params:mesh_tagging.email",
                    "api_key": "params:mesh_tagging.api_key",
                    "batch_name": "params:oa.data_processing.depth.mesh_tagging.batch_name",
                    "cmd": "params:oa.data_processing.depth.mesh_tagging.cmd",
                    "cmdargs": "params:oa.data_processing.depth.mesh_tagging.cmdargs",
                },
                outputs="oa.data_processing.depth.mesh_tagging.output",
            ),
            node(
                func=mesh_results_to_df,
                inputs={"mesh_results": "oa.data_processing.depth.mesh_tagging.output"},
                outputs="oa.data_processing.depth.mesh_tagging.intermediate",
            ),
        ],
        tags=["oa_depth_mesh_tagging"],
    )

    subfield_pipeline = pipeline(
        [
            node(
                func=process_subfield_data,
                inputs={
                    "data": "data_placeholder.raw",
                },
                outputs="data_placeholder.primary",
            )
        ],
    )

    subfield_pipelines = [
        pipeline(
            subfield_pipeline,
            inputs={
                "data_placeholder.raw": f"oa.data_collection.subfield.{subfield}.raw",
            },
            outputs={
                "data_placeholder.primary": f"oa.data_processing.subfield.{subfield}.primary",
            },
            tags=[
                f"oa.data_processing.subfield.{subfield}",
                "oa.data_processing.subfields",
            ],
        )
        for subfield in settings.DYNAMIC_PIPELINES_MAPPING["oa"]["subfields"]
    ]

    baseline_level0_pipeline = pipeline(
        [
            node(
                func=process_data_by_level,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.raw",
                    "level": "params:oa.data_collection.depth.levels.0",
                    "extra_mesh": "params:false_",
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
                    "extra_mesh": "params:false_",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.1.ptd.intermediate",
            ),
            node(
                func=concat_pq_ptd,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.1.ptd.intermediate",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.1.intermediate",
                tags="concat_pq_ptd"
            )
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
                outputs=["oa.data_processing.structural_biology.depth.reassigned.ct.intermediate", "oa.data_processing.structural_biology.depth.reassigned.other.intermediate"],
            )
        ],
        tags="oa.data_processing.depth.reassign_ct",
    )

    return (
        sum(additioal_mesh_levels)
        + sum(no_additional_mesh_levels)
        + mesh_tagging_pipeline
        + combine_levels_pipeline
        + sum(subfield_pipelines)
        + baseline_level0_pipeline
        + baseline_level1_pipeline
        + reassign_baseline_levels_pipeline
    )
