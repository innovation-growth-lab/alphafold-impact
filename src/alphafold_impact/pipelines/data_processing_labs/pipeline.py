"""
This is a boilerplate pipeline 'data_processing_labs'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import calculate_lab_determinants, combine_lab_results


def create_pipeline(  # pylint: disable=unused-argument&missing-function-docstring
    **kwargs,
) -> Pipeline:
    lab_determinants_pipeline = pipeline(
        [
            node(
                func=calculate_lab_determinants,
                inputs={
                    "dict_loader": "lab.data_collection.candidates.publications.intermediate"
                },
                outputs="lab.data_collection.candidates.scores.intermediate",
            )
        ],
        tags="data_processing_labs",
    )

    combine_lab_results_pipeline = pipeline(
        [
            node(
                func=combine_lab_results,
                inputs={
                    "dict_loader": "lab.data_collection.candidates.scores.intermediate"
                },
                outputs="lab.data_collection.candidates.scores.primary",
            )
        ],
        tags=["data_processing_labs", "combine_lab_results"],
    )

    return lab_determinants_pipeline + combine_lab_results_pipeline