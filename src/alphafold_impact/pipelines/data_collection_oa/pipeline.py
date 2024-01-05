from kedro.pipeline import Pipeline, node, pipeline
from .nodes import retrieve_oa_works_for_concepts_and_years


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=retrieve_oa_works_for_concepts_and_years,
                inputs=[
                    "params:test_concept_ids",
                    "params:test_publication_years",
                ],
                outputs="oa_raw_works_for_concepts_and_years",
            ),
        ]
    )
