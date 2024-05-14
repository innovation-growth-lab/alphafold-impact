"""
This is a boilerplate pipeline 'data_analysis_lab_productivity'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    get_sb_lab_outputs,
    compute_publication_production,
    get_event_study_outputs,
    get_event_study_strength,
    get_event_study_pdb_submissions,
)


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                get_sb_lab_outputs,
                inputs={
                    "data_loaders": "sb_lab.data_collection.publications.raw",
                    "mapping_df": "sb_lab.data_collection.candidates.map",
                },
                outputs="sb_lab.data_analysis.outputs.input",
                tags=["sb_lab_outputs"],
            ),
            node(
                compute_publication_production,
                inputs="sb_lab.data_analysis.outputs.input",
                outputs=[
                    "sb_lab.data_analysis.monthly_outputs",
                    "sb_lab.data_analysis.yearly_outputs",
                ],
                tags=["publication_production"],
            ),
            node(
                get_event_study_outputs,
                inputs={
                    "data": "sb_lab.data_analysis.outputs.input",
                    "level0": "analysis.descriptive.level0_data.processed",
                    "output_type": "params:analysis.output_type.publications",
                },
                outputs=[
                    "sb_lab.data_analysis.outputs.primary",
                    "sb_lab.data_analysis.outputs.event_study",
                ],
                tags=["event_study_outputs", "event_study"],
            ),
            node(
                get_event_study_strength,
                inputs={
                    "data": "sb_lab.data_analysis.outputs.primary",
                    "sc_data_af": "chains.complete_strong_links.id.primary",
                    "sc_data_ct": "chains.complete_strong_links.id.ct.primary",
                    "output_type": "params:analysis.output_type.publications",
                },
                outputs="sb_lab.data_analysis.outputs.event_study_strength",
                tags=["event_study_strength", "event_study"],
            ),
            node(
                get_event_study_pdb_submissions,
                inputs={
                    "data": "sb_lab.data_analysis.outputs.primary",
                    "pdb_data": "pdb.entries.intermediate",
                    "output_type": "params:analysis.output_type.publications",
                },
                outputs="sb_lab.data_analysis.outputs.event_study_pdb_submissions",
                tags=["event_study_pdb_submissions", "event_study"],
            ),
        ]
    )
