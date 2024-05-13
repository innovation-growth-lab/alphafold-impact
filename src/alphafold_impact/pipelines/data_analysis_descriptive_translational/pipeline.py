"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import load_input_data, merge_inputs, preprocess_sb_data, get_cc_counts


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                load_input_data,
                inputs={
                    "data": "oa.data_processing.depth.other.primary",
                    "source": "params:analysis.source.other",
                },
                outputs="analysis.other.level0_data",
                tags=["other_descriptive"]
            ),
            node(
                load_input_data,
                inputs={
                    "data": "oa.data_processing.depth.primary",
                    "source": "params:analysis.source.af",
                },
                outputs="analysis.af.level0_data",
                tags=["af_descriptive"]
            ),
            node(
                load_input_data,
                inputs={
                    "data": "oa.data_processing.depth.ct.primary",
                    "source": "params:analysis.source.ct",
                },
                outputs="analysis.ct.level0_data",
                tags=["ct_descriptive"]
            ),
            node(
                merge_inputs,
                inputs={
                    "alphafold_data": "analysis.af.level0_data",
                    "ct_data": "analysis.ct.level0_data",
                    "other_data": "analysis.other.level0_data",
                },
                outputs="analysis.descriptive.level0_data",
                tags="descriptive_merge"
            ),
            node(
                preprocess_sb_data,
                inputs="analysis.descriptive.level0_data",
                outputs="analysis.descriptive.level0_data.processed",
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
