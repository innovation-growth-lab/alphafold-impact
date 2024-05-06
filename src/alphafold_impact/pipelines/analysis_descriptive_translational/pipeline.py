"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import load_input_data, merge_inputs, get_cc_counts


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                load_input_data,
                inputs={
                    "data": "chains.complete_links.id.primary",
                    "source": "params:analysis.source.af",
                },
                outputs="af.level0_data",
            ),
            node(
                load_input_data,
                inputs={
                    "data": "oa.data_processing.structural_biology.depth.reassigned.ct.intermediate",
                    "source": "params:analysis.source.ct",
                },
                outputs="ct.level0_data",
            ),
            node(
                load_input_data,
                inputs={
                    "data": "oa.data_processing.structural_biology.depth.reassigned.other.intermediate",
                    "source": "params:analysis.source.other",
                },
                outputs="other.level0_data",
            ),
            node(
                merge_inputs,
                inputs={
                    "alphafold_data": "af.level0_data",
                    "ct_data": "ct.level0_data",
                    "other_data": "other.level0_data",
                },
                outputs="analysis.descriptive.level0_data",
            ),
            node(
                get_cc_counts,
                inputs={
                    "data": "analysis.descriptive.level0_data",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="analysis.descriptive.level0_data_with_cc_counts",
                tags=["cc_counts"]
            ),
        ],
        tags=["analysis_descriptive_translational"],
    )
