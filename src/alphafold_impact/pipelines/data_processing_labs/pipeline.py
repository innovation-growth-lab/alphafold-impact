"""
This is a boilerplate pipeline 'data_processing_labs'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import calculate_lab_determinants


def create_pipeline(  # pylint: disable=unused-argument&missing-function-docstring
    **kwargs,
) -> Pipeline:
    return pipeline(
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
