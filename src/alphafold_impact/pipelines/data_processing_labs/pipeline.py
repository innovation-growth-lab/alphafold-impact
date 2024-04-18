"""
This is a boilerplate pipeline 'data_processing_labs'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import calculate_lab_determinants, combine_lab_results, assign_lab_label


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

    assign_lab_label_pipeline = pipeline(
        [
            node(
                func=assign_lab_label,
                inputs={
                    "candidate_data": "lab.data_collection.candidates.scores.primary",
                    "ground_truth_data": "lab.data_collection.candidates.nih.intermediate",
                },
                outputs="lab.data_collection.assignment.primary",
            )
        ],
        tags=["data_processing_labs", "assign_lab_label"],
    )

    return lab_determinants_pipeline + combine_lab_results_pipeline + assign_lab_label_pipeline