"""Pipeline for data processing.

This pipeline processes NIH clinical trials data.
To run this pipeline, use the following command:

    $ kedro run --pipeline data_processing_nih_clinical_trials
"""
from kedro.pipeline import Pipeline, pipeline, node
from .nodes import trials_with_references


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=trials_with_references,
                inputs="nih.data_collection.clinical_trials.raw",
                outputs="nih.processing.clinical_trials_with_references.intermediate",
            ),
        ],
    )
