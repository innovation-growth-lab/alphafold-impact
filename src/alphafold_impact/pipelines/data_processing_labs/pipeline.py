"""
This is a boilerplate pipeline 'data_processing_labs'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    calculate_lab_determinants,
    combine_lab_results,
    assign_lab_label,
    get_discipline_map,
)


def create_pipeline(  # pylint: disable=unused-argument&missing-function-docstring
    **kwargs,
) -> Pipeline:
    candidate_mapping_pipeline = pipeline(
        [
            node(
                func=get_discipline_map,
                inputs={"dict_loader": "lab.data_collection.candidates.intermediate"},
                outputs="lab.data_collection.candidates.map.intermediate",
            )
        ],
        tags=["data_processing_labs", "labs_map"],
    )

    lab_determinants_pipeline = pipeline(
        [
            node(
                func=calculate_lab_determinants,
                inputs={
                    "dict_loader": "lab.data_collection.candidates.publications.intermediate",
                    "candidate_map": "lab.data_collection.candidates.map.intermediate",
                },
                outputs="lab.data_collection.candidates.scores.intermediate",
            )
        ],
        tags=["data_processing_labs", "calculate_determinants"],
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

    return (
        candidate_mapping_pipeline
        + lab_determinants_pipeline
        + combine_lab_results_pipeline
        + assign_lab_label_pipeline
    )
