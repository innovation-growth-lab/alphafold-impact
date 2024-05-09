"""
This is a boilerplate pipeline 'data_analysis_chains'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import process_af_data, process_baseline_data, assign_focal_label


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    return pipeline(
        [
            node(
                process_baseline_data,
                inputs={
                    "intent_data": "oa.data_processing.depth.level.0.baseline.primary",
                    "baseline_data": "oa.data_collection.subfield.structural_biology.depth.0.intermediate",
                    "seed_baseline_data": "oa.data_processing.subfield.structural_biology.primary",
                },
                outputs="baseline_candidates",
            ),
            node(
                process_af_data,
                inputs={"alphafold_data": "oa.data_processing.depth.all.primary"},
                outputs="alphafold_target",
            ),
            node(
                assign_focal_label,
                inputs={
                    "baseline_candidates": "baseline_candidates",
                    "alphafold_target": "alphafold_target",
                },
                outputs="chains.seed_technologies.intermediate",
            ),
        ],
        tags="data_analysis_chains",
    )
