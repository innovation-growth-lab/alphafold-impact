"""
This is a boilerplate pipeline 'data_processing_chains'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    filter_relevant_citation_links,
    get_papers_with_full_chain,
    get_papers_with_clinical_article_citations,
    get_chain_label_papers,
)


def create_pipeline(  # pylint: disable=unused-argument&missing-function-docstring
    **kwargs,
) -> Pipeline:
    filter_citation_links = pipeline(
        [
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.id",
                },
                outputs="chains.complete_links.id.primary",
                tags="chains.id",
            ),
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.pmid",
                },
                outputs="chains.complete_links.pmid.primary",
                tags="chains.pmid",
            ),
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.doi",
                },
                outputs="chains.complete_links.doi.primary",
                tags="chains.doi",
            ),
        ],
        tags="complete_chains",
    )

    strong_citation_links = pipeline(
        [
            node(
                get_papers_with_full_chain,
                inputs={
                    "chains": "chains.complete_links.id.primary",
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.id",
                },
                outputs="chains.complete_strong_links.id.primary",
                tags="chains.id",
            ),
            node(
                get_papers_with_full_chain,
                inputs={
                    "chains": "chains.complete_links.pmid.primary",
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.pmid",
                },
                outputs="chains.complete_strong_links.pmid.primary",
                tags="chains.pmid",
            ),
            node(
                get_papers_with_full_chain,
                inputs={
                    "chains": "chains.complete_links.doi.primary",
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.doi",
                },
                outputs="chains.complete_strong_links.doi.primary",
                tags="chains.doi",
            ),
        ],
        tags="strong_paper_chains",
    )

    clinical_article_citations = pipeline(
        [
            node(
                get_papers_with_clinical_article_citations,
                inputs={
                    "chains": "chains.complete_links.pmid.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                    "identifier": "params:chains.identifier.pmid",
                },
                outputs="chains.complete_strong_links_with_ca.pmid.primary",
                tags="chains.pmid",
            ),
            node(
                get_papers_with_clinical_article_citations,
                inputs={
                    "chains": "chains.complete_links.doi.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                    "identifier": "params:chains.identifier.doi",
                },
                outputs="chains.complete_strong_links_with_ca.doi.primary",
                tags="chains.doi",
            ),
        ],
        tags="clinical_chains",
    )

    label_dataset = pipeline(
        [
            node(
                get_chain_label_papers,
                inputs={
                    "chains": "chains.complete_links.id.primary",
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.id",
                },
                outputs="oa.chain_labels.id.primary",
                tags="label.id",
            )
        ],
        tags="label_dataset",
    )

    return (
        filter_citation_links
        + strong_citation_links
        + clinical_article_citations
        + label_dataset
    )