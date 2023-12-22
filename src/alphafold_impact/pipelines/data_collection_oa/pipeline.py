from kedro.pipeline import Pipeline, node, pipeline

from .nodes import oa_works_for_concepts_and_years


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=oa_works_for_concepts_and_years,
                inputs=[
                    "params:concept_ids",
                    "params:publication_years",
                ],
                outputs="oa_works",
                name="retrieve_oa_works_for_concepts_and_years_node",
            ),
        ]
    )
