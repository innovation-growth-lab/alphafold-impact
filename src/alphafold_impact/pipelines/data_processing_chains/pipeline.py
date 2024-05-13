"""
This is a boilerplate pipeline 'data_processing_chains'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    filter_relevant_citation_links,
    get_papers_with_strong_chain,
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
                tags=["chains.pmid", "test_af"],
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
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.id.primary",
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.id",
                },
                outputs="chains.complete_strong_links.id.primary",
                tags="chains.id",
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.pmid.primary",
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "identifier": "params:chains.identifier.pmid",
                },
                outputs="chains.complete_strong_links.pmid.primary",
                tags=["chains.pmid", "test_af"],
            ),
            node(
                get_papers_with_strong_chain,
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

    filter_citation_links_ct = pipeline(
        [
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.id",
                },
                outputs="chains.complete_links.id.ct.primary",
                tags="chains.id",
            ),
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.pmid",
                },
                outputs="chains.complete_links.pmid.ct.primary",
                tags="chains.pmid",
            ),
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.doi",
                },
                outputs="chains.complete_links.doi.ct.primary",
                tags="chains.doi",
            ),
        ],
        tags="complete_chains_ct.primary",
    )

    strong_citation_links_ct = pipeline(
        [
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.id.ct.primary",
                    "alphafold_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.id",
                },
                outputs="chains.complete_strong_links.id.ct.primary",
                tags=["chains.id", "test_ct"],
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.pmid.ct.primary",
                    "alphafold_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.pmid",
                },
                outputs="chains.complete_strong_links.pmid.ct.primary",
                tags="chains.pmid",
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.doi.ct.primary",
                    "alphafold_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.doi",
                },
                outputs="chains.complete_strong_links.doi.ct.primary",
                tags="chains.doi",
            ),
        ],
        tags="strong_paper_chains_ct.primary",
    )

    clinical_article_citations_ct = pipeline(
        [
            node(
                get_papers_with_clinical_article_citations,
                inputs={
                    "chains": "chains.complete_links.pmid.ct.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                    "identifier": "params:chains.identifier.pmid",
                },
                outputs="chains.complete_strong_links_with_ca.pmid.ct.primary",
                tags="chains.pmid",
            ),
            node(
                get_papers_with_clinical_article_citations,
                inputs={
                    "chains": "chains.complete_links.doi.ct.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                    "identifier": "params:chains.identifier.doi",
                },
                outputs="chains.complete_strong_links_with_ca.doi.ct.primary",
                tags="chains.doi",
            ),
        ],
        tags="clinical_chains_ct.primary",
    )

    label_dataset_ct = pipeline(
        [
            node(
                get_chain_label_papers,
                inputs={
                    "chains": "chains.complete_links.id.ct.primary",
                    "alphafold_data": "oa.data_processing.depth.ct.primary",
                    "identifier": "params:chains.identifier.id",
                },
                outputs="oa.chain_labels.id.ct.primary",
                tags="label.id",
            )
        ],
        tags="label_dataset_ct.primary",
    )

    filter_citation_links_other = pipeline(
        [
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.id",
                    "num_levels": "params:chains.num_levels.other"
                },
                outputs="chains.complete_links.id.other.primary",
                tags="chains.id",
            ),
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.pmid",
                    "num_levels": "params:chains.num_levels.other"
                },
                outputs="chains.complete_links.pmid.other.primary",
                tags="chains.pmid",
            ),
            node(
                filter_relevant_citation_links,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.doi",
                    "num_levels": "params:chains.num_levels.other"
                },
                outputs="chains.complete_links.doi.other.primary",
                tags="chains.doi",
            ),
        ],
        tags=["complete_chains_other.primary", "chains_other"],
    )

    strong_citation_links_other = pipeline(
        [
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.id.other.primary",
                    "alphafold_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.id"
                },
                outputs="chains.complete_strong_links.id.other.primary",
                tags=["chains.id", "test_other"],
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.pmid.other.primary",
                    "alphafold_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.pmid",
                },
                outputs="chains.complete_strong_links.pmid.other.primary",
                tags="chains.pmid",
            ),
            node(
                get_papers_with_strong_chain,
                inputs={
                    "chains": "chains.complete_links.doi.other.primary",
                    "alphafold_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.doi",
                },
                outputs="chains.complete_strong_links.doi.other.primary",
                tags="chains.doi",
            ),
        ],
        tags=["strong_paper_chains_other.primary", "chains_other"],
    )

    clinical_article_citations_other = pipeline(
        [
            node(
                get_papers_with_clinical_article_citations,
                inputs={
                    "chains": "chains.complete_links.pmid.other.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                    "identifier": "params:chains.identifier.pmid",
                },
                outputs="chains.complete_strong_links_with_ca.pmid.other.primary",
                tags="chains.pmid",
            ),
            node(
                get_papers_with_clinical_article_citations,
                inputs={
                    "chains": "chains.complete_links.doi.other.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                    "identifier": "params:chains.identifier.doi",
                },
                outputs="chains.complete_strong_links_with_ca.doi.other.primary",
                tags="chains.doi",
            ),
        ],
        tags=["clinical_chains_other.primary", "chains_other"],
    )

    label_dataset_other = pipeline(
        [
            node(
                get_chain_label_papers,
                inputs={
                    "chains": "chains.complete_links.id.other.primary",
                    "alphafold_data": "oa.data_processing.depth.other.primary",
                    "identifier": "params:chains.identifier.id",
                },
                outputs="oa.chain_labels.id.other.primary",
                tags="label.id",
            )
        ],
        tags=["label_dataset_other.primary", "chains_other"],
    )

    return (
        filter_citation_links
        + strong_citation_links
        + clinical_article_citations
        + label_dataset
        + filter_citation_links_ct
        + strong_citation_links_ct
        + clinical_article_citations_ct
        + label_dataset_ct
        + filter_citation_links_other
        + strong_citation_links_other
        + clinical_article_citations_other
        + label_dataset_other
    )
