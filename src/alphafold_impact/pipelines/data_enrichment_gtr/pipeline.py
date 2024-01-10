"""
Pipeline for collecting OpenAlex publications from that correspond to a
given Gateway to Research publication.

To run this pipeline, use the following command:

    $ kedro run --pipeline data_enrichment_gtr
"""

from kedro.pipeline import Pipeline, pipeline, node
from alphafold_impact.pipelines.data_collection_oa.nodes import collect_papers
from .nodes import (
    preprocess_publication_doi,
    create_list_doi_inputs,
    create_dictionary_doi_to_oa,
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    return pipeline(
        [
            node(
                func=preprocess_publication_doi,
                inputs="gtr_intermediate_publications",
                outputs="gtr_preprocessed_publications",
            ),
            node(
                func=create_list_doi_inputs,
                inputs="gtr_preprocessed_publications",
                outputs="gtr_doi_list",
            ),
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:api.mailto",
                    "perpage": "params:api.perpage",
                    "work_ids": "gtr_doi_list",
                    "filter_by": "params:oa_gtr.filter",
                    "group_work_ids": "params:oa.group_work_ids",
                    "eager_loading": "params:oa.eager_loading",
                },
                outputs="gtr_dict_matches",
            ),
            node(
                func=create_dictionary_doi_to_oa,
                inputs="gtr_dict_matches",
                outputs="gtr_intermediate_publications_oa_dict",
            ),
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:api.mailto",
                    "perpage": "params:api.perpage",
                    "work_ids": "gtr_intermediate_publications_oa_dict",
                    "filter_by": "params:direction.cited_by", # 
                },
                outputs="gtr_intermediate_publication_citations",
            ),
        ]
    )
