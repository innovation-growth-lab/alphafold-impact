"""
This is a boilerplate pipeline 'data_processing_oa'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, node, pipeline
from alphafold_impact import settings
from ...utils.nih_mesh_tagging import (  # pylint: disable=E0402
    generate_text_for_mesh_tagging,
    skr_web_python_api_generic_batch,
    mesh_results_to_df,
)
from .nodes import (
    combine_depth_strength_level_0, 
    combine_depth_strength_other_levels,
    process_subfield_data
)


def create_pipeline(**kwargs) -> Pipeline: # pylint: disable=unused-argument&missing-function-docstring

    merge_zeroth_level = pipeline(
        [
            node(
                combine_depth_strength_level_0,
                inputs={
                    "depth_data": "oa.data_collection.depth.no_mesh.0.intermediate",
                    "strength_data": "s2.data_collection.strength.level.0",
                },
                outputs="oa.data_processing.depth.level.0.primary",
                name="merge_zeroth_level",
            )
        ],
        tags="oa.data_processing.depth.add_strength",
    )

    merge_pipelines = [
        pipeline(
            [
                node(
                    combine_depth_strength_other_levels,
                    inputs={
                        "previous_depth_data": f"oa.data_collection.depth.no_mesh.{level-1}.intermediate",
                        "depth_data": f"oa.data_collection.depth.no_mesh.{level}.intermediate",
                        "strength_data": f"s2.data_collection.strength.level.{level}",
                    },
                    outputs=f"oa.data_processing.depth.level.{level}.primary",
                    name=f"merge_level_{level}",
                )
            ],
            tags="oa.data_processing.depth.add_strength",
        )
        for level in settings.DYNAMIC_PIPELINES_MAPPING["depth_levels"][1:]
    ]

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
                    "data": "oa.data_collection.subfield.placeholder.raw",
                },
                outputs="oa.data_processing.subfield.placeholder.primary",
            )
        ],
    )

    subfield_pipelines = [
        pipeline(
            subfield_pipeline,
            inputs={
                "oa.data_collection.subfield.placeholder.raw": f"oa.data_collection.subfield.{subfield}.raw",
            },
            outputs={
                "oa.data_processing.subfield.placeholder.primary": f"oa.data_processing.subfield.{subfield}.primary",
            },
            tags=[
                f"oa.data_processing.subfield.{subfield}",
                "oa.data_processing.subfields",
            ]
        )
        for subfield in settings.DYNAMIC_PIPELINES_MAPPING["oa"]["subfields"]
    ]

    return mesh_tagging_pipeline + merge_zeroth_level + sum(merge_pipelines) + sum(subfield_pipelines)
