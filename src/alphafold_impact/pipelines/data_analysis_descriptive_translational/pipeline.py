"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    load_input_data,
    load_input_applied_data,
    merge_inputs,
    preprocess_sb_data,
    get_cc_papers,
    get_cc_moderators,
    get_patent_papers,
    get_patent_moderators,
    get_patent_classifications,
    create_tcc_sb_papers,
    create_publications_data,
    merge_individual_data
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    level0_pipeline = pipeline(
        [
            node(
                load_input_data,
                inputs={
                    "data": "oa.data_processing.depth.all.primary",
                    "source": "params:analysis.source.af",
                },
                outputs="analysis.af.level0_data",
                tags=["af_descriptive"],
            ),
            node(
                load_input_data,
                inputs={
                    "data": "oa.data_processing.depth.ct.primary",
                    "source": "params:analysis.source.ct",
                },
                outputs="analysis.ct.level0_data",
                tags=["ct_descriptive"],
            ),
            node(
                load_input_data,
                inputs={
                    "data": "oa.data_processing.depth.other.primary",
                    "source": "params:analysis.source.other",
                },
                outputs="analysis.other.level0_data",
                tags=["other_descriptive"],
            ),
            node(
                merge_inputs,
                inputs={
                    "alphafold_data": "analysis.af.level0_data",
                    "ct_data": "analysis.ct.level0_data",
                    "other_data": "analysis.other.level0_data",
                },
                outputs="analysis.descriptive.level0_data",
                tags="descriptive_merge",
            ),
            node(
                preprocess_sb_data,
                inputs="analysis.descriptive.level0_data",
                outputs="analysis.descriptive.level0_data.processed",
                tags="descriptive_preprocess",
            ),
            node(
                get_cc_papers,
                inputs={
                    "data": "analysis.descriptive.level0_data.processed",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="analysis.descriptive.level0_data_with_cc",
                tags=["cc_counts"],
            ),
            node(
                get_cc_moderators,
                inputs={
                    "cc_data": "analysis.descriptive.level0_data_with_cc",
                    "sc_data_af": "chains.complete_strong_links.id.primary",
                    "sc_data_ct": "chains.complete_strong_links.id.ct.primary",
                    "pdb_data": "pdb.entries.intermediate",
                },
                outputs="analysis.descriptive.level0_data_with_moderators_cc",
                tags=["cc_moderators", "moderators"],
            ),
            node(
                get_patent_papers,
                inputs={
                    "data": "analysis.descriptive.level0_data.processed",
                    "patent_data": "lens.data_processing.primary",
                },
                outputs="analysis.descriptive.level0_data_with_patents",
                tags=["patent_counts"],
            ),
            node(
                get_patent_moderators,
                inputs={
                    "pc_data": "analysis.descriptive.level0_data_with_patents",
                    "sc_data_af": "chains.complete_strong_links.id.primary",
                    "sc_data_ct": "chains.complete_strong_links.id.ct.primary",
                    "pdb_data": "pdb.entries.intermediate",
                },
                outputs="analysis.descriptive.level0_data_with_moderators_pc",
                tags=["pc_moderators", "moderators"],
            ),
        ],
        tags=["analysis_descriptive_translational"],
    )

    applied_pipeline = pipeline(
        [
            node(
                load_input_applied_data,
                inputs={
                    "data": "oa.data_processing.depth.all.primary",
                    "source": "params:analysis.source.af",
                },
                outputs="analysis.af.applied_data",
                tags=["af_descriptive", "applied_strength_descriptive"],
            ),
            node(
                load_input_applied_data,
                inputs={
                    "data": "oa.data_processing.depth.ct.primary",
                    "source": "params:analysis.source.ct",
                },
                outputs="analysis.ct.applied_data",
                tags=["ct_descriptive", "applied_strength_descriptive"],
            ),
            node(
                load_input_applied_data,
                inputs={
                    "data": "oa.data_processing.depth.other.primary",
                    "source": "params:analysis.source.other",
                },
                outputs="analysis.other.applied_data",
                tags=["other_descriptive", "applied_strength_descriptive"],
            ),
            node(
                merge_inputs,
                inputs={
                    "alphafold_data": "analysis.af.applied_data",
                    "ct_data": "analysis.ct.applied_data",
                    "other_data": "analysis.other.applied_data",
                },
                outputs="analysis.descriptive.applied_data",
                tags=["descriptive_merge", "applied_strength_descriptive"],
            ),
            node(
                get_cc_papers,
                inputs={
                    "data": "analysis.descriptive.applied_data",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="analysis.descriptive.applied_data_with_cc",
                tags=["cc_counts", "cc_applied"],
            ),
            node(
                get_cc_moderators,
                inputs={
                    "cc_data": "analysis.descriptive.applied_data_with_cc",
                    "sc_data_af": "chains.complete_strong_links.id.primary",
                    "sc_data_ct": "chains.complete_strong_links.id.ct.primary",
                    "pdb_data": "pdb.entries.intermediate",
                },
                outputs="analysis.descriptive.applied_data_with_moderators_cc",
                tags=["cc_moderators", "moderators", "cc_applied"],
            ),
            node(
                get_patent_papers,
                inputs={
                    "data": "analysis.descriptive.applied_data",
                    "patent_data": "lens.data_processing.primary",
                },
                outputs="analysis.descriptive.applied_data_with_patents",
                tags=["patent_counts", "applied_strength_descriptive", "cc_applied"],
            ),
            node(
                get_patent_moderators,
                inputs={
                    "pc_data": "analysis.descriptive.applied_data_with_patents",
                    "sc_data_af": "chains.complete_strong_links.id.primary",
                    "sc_data_ct": "chains.complete_strong_links.id.ct.primary",
                    "pdb_data": "pdb.entries.intermediate",
                },
                outputs="analysis.descriptive.applied_data_with_moderators_pc",
                tags=["pc_moderators", "moderators", "cc_applied"],
            ),
        ],
        tags=["analysis_descriptive_translational_applied"],
    )

    patent_cpc_pipeline = pipeline(
        [
            node(
                get_patent_classifications,
                inputs={
                    "patent_data": "analysis.descriptive.level0_data_with_patents",
                    "cpc_codes": "cpc.codes",
                },
                outputs="analysis.descriptive.level0_data_with_patents_cpc",
                tags=["patent_cpc"],
            ),
        ],
        tags=["analysis_descriptive_translational_patent_cpc"],
    )

    tcc_sb_papers = pipeline(
        [
            node(
                create_tcc_sb_papers,
                inputs={
                    "data": "analysis.descriptive.level0_data_with_cc",
                },
                outputs=["analysis.descriptive.tcc_sb_papers", "reporting.tcc_sb_papers"],
                tags=["tcc_sb_papers"],
            ),
        ],
        tags=["analysis_descriptive_translational_tcc_sb_papers"],
    )

    create_publications_data_pipeline = pipeline(
        [
            node(
                create_publications_data,
                inputs={
                    "data": "oa.chain_labels.id.primary",
                    "source": "params:analysis.source.af",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "patents_data": "lens.data_processing.primary",
                    "grants_data": "oa.data_processing.depth.grants.primary",
                    "pdb_submissions": "pdb.entries.intermediate",
                },
                outputs="analysis.descriptive.data.af",
                tags=["af_publications_data"],
            ),
            node(
                create_publications_data,
                inputs={
                    "data": "oa.chain_labels.id.ct.primary",
                    "source": "params:analysis.source.ct",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "patents_data": "lens.data_processing.primary",
                    "grants_data": "oa.data_processing.depth.grants.primary",
                    "pdb_submissions": "pdb.entries.intermediate",
                },
                outputs="analysis.descriptive.data.ct",
                tags=["ct_publications_data"],
            ),
            node(
                create_publications_data,
                inputs={
                    "data": "oa.data_processing.structural_biology.depth.other.intermediate",
                    "source": "params:analysis.source.other",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "patents_data": "lens.data_processing.primary",
                    "grants_data": "oa.data_processing.depth.grants.primary",
                    "pdb_submissions": "pdb.entries.intermediate",
                },
                outputs="analysis.descriptive.data.other",
                tags=["other_publications_data"],
            ),
            node(
                merge_individual_data,
                inputs={
                    "data_af": "analysis.descriptive.data.af",
                    "data_ct": "analysis.descriptive.data.ct",
                    "data_other": "analysis.descriptive.data.other",
                },
                outputs="analysis.descriptive.data.outputs",
                tags=["merge_publications_data"],
            ),
        ],
        tags=["create_individual_data"],
    )

    return level0_pipeline + applied_pipeline + patent_cpc_pipeline + tcc_sb_papers + create_publications_data_pipeline
