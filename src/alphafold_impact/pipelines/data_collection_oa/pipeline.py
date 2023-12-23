from kedro.pipeline import Pipeline, node, pipeline

from .nodes import retrieve_oa_works_for_concepts_and_years, save_oa_works_to_s3


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=retrieve_oa_works_for_concepts_and_years,
                inputs=[
                    "params:concept_ids",
                    "params:publication_years",
                ],
                outputs="oa_works_for_concepts_list",
            ),
            node(
                func=save_oa_works_to_s3,
                inputs="oa_works_for_concepts_list",
                outputs="oa_works_for_concepts_json",
            ),
        ]
    )
