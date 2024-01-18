"""
Pipeline for collecting OpenAlex publications that match the doi values
for the Gateway to Research publications. The pipeline then collects all
papers that these publications cite.

To run this pipeline, use the following command:

    $ kedro run --pipeline data_enrichment_gtr
"""

from kedro.pipeline import Pipeline, pipeline, node
from alphafold_impact.pipelines.data_collection_oa.nodes import collect_papers
from .nodes import (
    preprocess_publication_doi,
    create_list_doi_inputs,
    test_node,
    load_referenced_work_ids,
)

LABEL = "gtr"

def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    return pipeline(
        [
            node(
                func=preprocess_publication_doi,
                inputs="gtr.collection.outcomes/publications.intermediate",
                outputs=f"{LABEL}.enrichment.publications.preproc",
                tags="gtr.publications",
            ),
            node(
                func=create_list_doi_inputs,
                inputs=f"{LABEL}.enrichment.publications.preproc",
                outputs=f"{LABEL}.enrichment.publications.doi.list",
                tags="gtr.publications",
            ),
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:oa.api.mailto",
                    "perpage": "params:oa.api.perpage",
                    "work_ids": f"{LABEL}.enrichment.publications.doi.list",
                    "filter_by": "params:oa.gtr.filter_doi",
                    "group_work_ids": "params:oa.group_work_ids",
                    "timestamp_keys": "params:oa.timestamp_keys",
                    "parallelise": "params:oa.parallelise",
                },
                outputs=f"{LABEL}.enrichment.publications.oa.tmp",
                tags="gtr.publications",
            ),
            node(
                func=load_referenced_work_ids,
                inputs=f"{LABEL}.enrichment.publications.oa.tmp",
                outputs=[
                    f"{LABEL}.enrichment.publications.oa.list",
                    f"{LABEL}.enrichment.publications.oa.dict",
                ],
                tags="gtr.citations",
            ),
            # node(
            #     func=create_dictionary_doi_to_oa,
            #     inputs="gtr.publications.oa.list",
            #     outputs="gtr.publications.dict",
            # ),
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:oa.api.mailto",
                    "perpage": "params:oa.api.perpage",
                    "work_ids": f"{LABEL}.enrichment.publications.oa.list",
                    "filter_by": "params:oa.gtr.filter_oa",
                    "group_work_ids": "params:oa.group_work_ids",
                    "timestamp_keys": "params:oa.timestamp_keys",
                    "parallelise": "params:oa.parallelise",
                },
                outputs=f"{LABEL}.enrichment.publications.oa.citations",
                tags="gtr.citations",
            ),
        ],
    )
