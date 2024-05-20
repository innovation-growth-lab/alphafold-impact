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
                tags=["af_descriptive"],
            ),
            node(
                load_input_applied_data,
                inputs={
                    "data": "oa.data_processing.depth.ct.primary",
                    "source": "params:analysis.source.ct",
                },
                outputs="analysis.ct.applied_data",
                tags=["ct_descriptive"],
            ),
            node(
                load_input_applied_data,
                inputs={
                    "data": "oa.data_processing.depth.other.primary",
                    "source": "params:analysis.source.other",
                },
                outputs="analysis.other.applied_data",
                tags=["other_descriptive"],
            ),
            node(
                merge_inputs,
                inputs={
                    "alphafold_data": "analysis.af.applied_data",
                    "ct_data": "analysis.ct.applied_data",
                    "other_data": "analysis.other.applied_data",
                },
                outputs="analysis.descriptive.applied_data",
                tags="descriptive_merge",
            ),
            node(
                get_cc_papers,
                inputs={
                    "data": "analysis.descriptive.applied_data",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="analysis.descriptive.applied_data_with_cc",
                tags=["cc_counts"],
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
                tags=["cc_moderators", "moderators"],
            ),
            node(
                get_patent_papers,
                inputs={
                    "data": "analysis.descriptive.applied_data",
                    "patent_data": "lens.data_processing.primary",
                },
                outputs="analysis.descriptive.applied_data_with_patents",
                tags=["patent_counts"],
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
                tags=["pc_moderators", "moderators"],
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

    return level0_pipeline + applied_pipeline + patent_cpc_pipeline
