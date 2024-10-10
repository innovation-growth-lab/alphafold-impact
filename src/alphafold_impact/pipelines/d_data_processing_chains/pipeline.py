"""
This is a boilerplate pipeline 'data_processing_chains'
generated using Kedro 0.19.1
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
    af_chain_pipeline = pipeline(
        [
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.af",
                },
                outputs="chains.complete_links.id.primary",
                name="clean_citation_chains_af"
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.id.primary",
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.af",
                },
                outputs="chains.complete_strong_links.id.primary",
                name="filter_strong_chains_af"
            ),
            node(
                get_chain_label_papers,
                inputs={
                    "chains": "chains.complete_links.id.primary",
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.af",
                },
                outputs="oa.chain_labels.id.primary",
                name="get_chain_labels_af"
            ),
        ],
        tags=["data_processing_chains"],
    )

    ct_chain_pipeline = pipeline(
        [
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.ct",
                },
                outputs="chains.complete_links.id.ct.primary",
                name="clean_citation_chains_ct"
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.id.ct.primary",
                    "alphafold_data": "oa.data_processing.depth.ct.primary",
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
                    "alphafold_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.ct",
                },
                outputs="oa.chain_labels.id.ct.primary",
                name="get_chain_labels_ct"
            ),
        ],
        tags=["data_processing_chains"],
    )

    other_chain_pipeline = pipeline(
        [
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.other",
                },
                outputs="chains.complete_links.id.other.primary",
                name="clean_citation_chains_other"
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.id.other.primary",
                    "alphafold_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.other",
                },
                outputs="chains.complete_strong_links.id.other.primary",
                name="filter_strong_chains_other"
            ),
            node(
                get_chain_label_papers,
                inputs={
                    "chains": "chains.complete_links.id.other.primary",
                    "alphafold_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.other",
                },
                outputs="oa.chain_labels.id.other.primary",
                name="get_chain_labels_other"
            ),
        ],
        tags=["data_processing_chains"],
    )

    return af_chain_pipeline + ct_chain_pipeline + other_chain_pipeline
