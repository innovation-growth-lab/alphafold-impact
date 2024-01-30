"""Pipeline for data collection.

This pipeline fetches NIH clinical trials data.
To run this pipeline, use the following command:

    $ kedro run --pipeline data_collection_gtr
"""
from kedro.pipeline import Pipeline, pipeline, node
from .nodes import collect_and_save_nih_clinical_trials


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=collect_and_save_nih_clinical_trials,
                inputs={
                    "nih_ct_download_url": "params:download_url",
                    "nih_ct_zip_save_path": "params:zip_save_path",
                    "nih_ct_extract_path": "params:extract_path",
                },
                outputs="raw",
            ),
        ],
        namespace="nih.data_collection.clinical_trials",
    )
