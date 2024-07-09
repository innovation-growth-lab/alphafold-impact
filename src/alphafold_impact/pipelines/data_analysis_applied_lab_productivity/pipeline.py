"""
This is a boilerplate pipeline 'data_analysis_lab_productivity'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    get_applied_lab_outputs,
    preprocess_for_event_study,
    get_applied_lab_staggered_outputs,
)
from ..data_analysis_lab_productivity.nodes import (  # pylint: disable=E0402
    compute_publication_production,
    # preprocess_for_event_study,
    get_event_study_outputs,
    get_event_study_strength,
    get_event_study_pdb_submissions,
    get_event_study_predictive_outputs,
    get_event_study_cc,
    get_event_study_pc,
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    basic_pipeline = pipeline(
        [
            node(
                get_applied_lab_outputs,
                inputs={
                    "data_loaders": "other_lab.data_collection.publications.raw",
                    "mapping_df": "other_lab.data_collection.candidates.map",
                },
                outputs="applied_lab.data_analysis.outputs.input",
                tags=["applied_lab_outputs", "event_study", "applied_cc", "event_study_strength_applied", "event_study_pc"],
            ),
            node(
                compute_publication_production,
                inputs="applied_lab.data_analysis.outputs.input",
                outputs=[
                    "applied_lab.data_analysis.monthly_outputs",
                    "applied_lab.data_analysis.yearly_outputs",
                ],
                tags=["publication_production"],
            ),
            node(
                preprocess_for_event_study,
                inputs={
                    "data": "applied_lab.data_analysis.outputs.input",
                    "applied_levels": "analysis.descriptive.applied_data",
                },
                outputs="applied_lab.data_analysis.outputs.primary",
                tags=["preprocess_for_event_study", "event_study", "event_study_pc", "applied_cc", "event_study_strength_applied"],
            ),
            node(
                get_event_study_outputs,
                inputs={
                    "data": "applied_lab.data_analysis.outputs.primary",
                },
                outputs=[
                    "applied_lab.data_analysis.outputs.counts.event_study",
                    "applied_lab.data_analysis.outputs.citations.event_study",
                ],
                tags=["event_study_outputs", "event_study"],
            ),
            node(
                get_event_study_strength,
                inputs={
                    "data": "applied_lab.data_analysis.outputs.primary",
                    "sc_data_af": "chains.complete_strong_links.id.primary",
                    "sc_data_ct": "chains.complete_strong_links.id.ct.primary",
                },
                outputs=[
                    "applied_lab.data_analysis.outputs.counts.event_study_strength",
                    "applied_lab.data_analysis.outputs.citations.event_study_strength",
                ],
                tags=["event_study_strength_applied", "event_study"],
            ),
            node(
                get_event_study_pdb_submissions,
                inputs={
                    "data": "applied_lab.data_analysis.outputs.primary",
                    "pdb_data": "pdb.entries.intermediate",
                },
                outputs=[
                    "applied_lab.data_analysis.outputs.counts.event_study_pdb_submissions",
                    "applied_lab.data_analysis.outputs.citations.event_study_pdb_submissions",
                ],
                tags=["event_study_pdb_submissions", "event_study"],
            ),
            node(
                get_event_study_predictive_outputs,
                inputs={
                    "data": "applied_lab.data_analysis.outputs.primary",
                },
                outputs=[
                    "applied_lab.data_analysis.outputs.counts.event_study_predictive",
                    "applied_lab.data_analysis.outputs.citations.event_study_predictive",
                ],
                tags=["event_study_predictive", "event_study"],
            ),
            node(
                get_event_study_cc,
                inputs={
                    "data": "applied_lab.data_analysis.outputs.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs=["applied_lab.data_analysis.outputs.papers_with_ccs",
                    "applied_lab.data_analysis.outputs.counts.event_study_cc",
                    "applied_lab.data_analysis.outputs.citations.event_study_cc"],
                
                tags=["event_study_cc", "event_study", "applied_cc"],
            ),
            node(
                get_event_study_pc,
                inputs={
                    "data": "applied_lab.data_analysis.outputs.primary",
                    "patent_data": "lens.data_processing.primary",
                },
                outputs=[
                    "applied_lab.data_analysis.outputs.papers_with_pcs",
                    "applied_lab.data_analysis.outputs.counts.event_study_pc",
                    "applied_lab.data_analysis.outputs.citations.event_study_pc",
                ],
                tags=["event_study_pc", "event_study", "pc"],
            ),
        ],
        tags="applied_lab_productivity" 
    )

    staggered_pipeline = pipeline(
        [
            node(
                get_applied_lab_staggered_outputs,
                inputs={
                    "data_loaders": "other_lab.data_collection.publications.raw",
                    "mapping_df": "other_lab.data_collection.candidates.map",
                    "sb_mapping_df": "sb_lab.data_collection.candidates.map",
                    "levels": "analysis.descriptive.applied_data",
                    "pdb_submissions": "pdb.entries.intermediate",
                    "strength_es": "applied_lab.data_analysis.outputs.counts.event_study_strength",
                    "patents_data": "lens.data_processing.primary",
                    "clinical_citations": "applied_lab.data_analysis.outputs.papers_with_ccs",
                    "grants_data": "oa.data_processing.depth.grants.primary",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "institutional_data": "other_lab.data_collection.institution_info.primary"
                },
                outputs= [
                    "applied_lab.data_analysis.staggered.outputs.primary",
                    "applied_lab.data_analysis.staggered.outputs.quarterly.primary",
                    "applied_lab.data_analysis.staggered.outputs.collapsed.primary",
                ],
                tags=["applied_lab_outputs"],
            ),
        ],
        tags=["staggered_applied", "staggered"],
    )

    return basic_pipeline + staggered_pipeline
