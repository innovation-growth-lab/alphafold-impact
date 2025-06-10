"""
This is the final pipeline for processing citation chains and generating the ultimate outputs.

PIPELINE FLOW:
    STEP 7: This pipeline (d_data_processing_chains) performs final citation chain analysis
    FROM: ← c_data_collection_s2 (provides data with citation intent from Semantic Scholar)
    FINAL OUTPUT: Chain labels and analysis for AlphaFold, counterfactual, and other papers

    This is the FINAL step in the pipeline sequence:
    a_data_collection_oa → b__data_processing_oa → b_data_processing_baselines →
    a_data_collection_oa → b__data_processing_oa → c_data_collection_s2 → d_data_processing_chains
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    filter_relevant_citation_links,
    get_papers_with_strong_chain,
    get_chain_label_papers,
)


def create_pipeline(  # pylint: disable=unused-argument&missing-function-docstring
    **kwargs,
) -> Pipeline:
    # STEP 7A: Analyse AlphaFold citation chains
    # ← FROM: c_data_collection_s2 (oa.data_processing.depth.af.primary)
    # FINAL OUTPUT: Chain analysis for AlphaFold papers
    af_chain_pipeline = pipeline(
        [
            node(
                filter_relevant_citation_links,
                inputs={
                    "depth_data": "oa.data_processing.depth.af.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.af",
                },
                outputs="chains.complete_links.af.primary",
                name="clean_citation_chains_af",
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.af.primary",
                    "depth_data": "oa.data_processing.depth.af.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.af",
                },
                outputs="chains.complete_strong_links.af.primary",
                name="filter_strong_chains_af",
            ),
            node(
                get_chain_label_papers,
                inputs={
                    "chains": "chains.complete_links.af.primary",
                    "depth_data": "oa.data_processing.depth.af.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.af",
                },
                outputs="oa.chain_labels.af.primary",
                name="get_chain_labels_af",
            ),
        ],
        tags=["data_processing_af_chains"],
    )

    # STEP 7B: Analyse counterfactual citation chains
    # ← FROM: c_data_collection_s2 (oa.data_processing.depth.ct.primary)
    # FINAL OUTPUT: Chain analysis for counterfactual papers
    ct_chain_pipeline = pipeline(
        [
            node(
                filter_relevant_citation_links,
                inputs={
                    "depth_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.ct",
                },
                outputs="chains.complete_links.id.ct.primary",
                name="clean_citation_chains_ct",
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.id.ct.primary",
                    "depth_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.ct",
                },
                outputs="chains.complete_strong_links.id.ct.primary",
                name="filter_strong_chains_ct",
            ),
            node(
                get_chain_label_papers,
                inputs={
                    "chains": "chains.complete_links.id.ct.primary",
                    "depth_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.ct",
                },
                outputs="oa.chain_labels.id.ct.primary",
                name="get_chain_labels_ct",
            ),
        ],
        tags=["data_processing_ct_chains"],
    )

    # STEP 7C: Analyse other structural biology citation chains
    # ← FROM: c_data_collection_s2 (oa.data_processing.depth.other.primary)
    # FINAL OUTPUT: Chain analysis for other structural biology papers
    other_chain_pipeline = pipeline(
        [
            node(
                filter_relevant_citation_links,
                inputs={
                    "depth_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.other",
                },
                outputs="chains.complete_links.id.other.primary",
                name="clean_citation_chains_other",
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.id.other.primary",
                    "depth_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.other",
                },
                outputs="chains.complete_strong_links.id.other.primary",
                name="filter_strong_chains_other",
            ),
            node(
                get_chain_label_papers,
                inputs={
                    "chains": "chains.complete_links.id.other.primary",
                    "depth_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.other",
                },
                outputs="oa.chain_labels.id.other.primary",
                name="get_chain_labels_other",
            ),
        ],
        tags=["data_processing_other_chains"],
    )

    return af_chain_pipeline + ct_chain_pipeline + other_chain_pipeline
