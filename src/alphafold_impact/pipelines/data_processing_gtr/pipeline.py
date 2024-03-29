# """
# Pipeline for collecting OpenAlex publications that match the doi values
# for the Gateway to Research publications. The pipeline then collects all
# papers that these publications cite.

# To run this pipeline, use the following command:

#     $ kedro run --pipeline data_enrichment_gtr
# """

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    load_gtr_data,
    process_institutions,
    process_publications,
    process_alphafold_citations
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    return pipeline(
        [
            node(
                func=load_gtr_data,
                inputs=["works"],
                outputs=["authors.raw", "publications.raw"],
                tags=["gtr.processing", "gtr.institutions", "gtr.publications"],
            ),
            node(
                func=process_institutions,
                inputs=["authors.raw"],
                outputs="institutions",
                tags=["gtr.processing", "gtr.institutions"],
            ),
            node(
                func=process_publications,
                inputs=["publications.raw"],
                outputs="publications",
                tags=["gtr.processing", "gtr.publications"],
            ),
            node(
                func=process_alphafold_citations,
                inputs=["dict", "cites.raw"],
                outputs="citations.alphafold",
                tags=["gtr.processing", "gtr.citations"],
            ),
        ],
        namespace="oa.data_processing.gtr"
    )
