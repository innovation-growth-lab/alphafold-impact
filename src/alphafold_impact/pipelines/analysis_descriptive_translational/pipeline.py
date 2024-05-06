"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import load_input_data


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                load_input_data,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "ct_data": "oa.data_processing.structural_biology.depth.reassigned.ct.intermediate",
                    "other_data": "oa.data_processing.structural_biology.depth.reassigned.other.intermediate",
                },
                outputs=["analysis.descriptive.level0_data"],
            ),
        ],
        tags=["analysis_descriptive_translational"],
    )
