"""Pipeline for PubMed and iCite data collection.

This pipeline provides two separate data collection workflows:
1. ClinicalTrials.gov data collection using the API v2 endpoints
2. NIH iCite data collection from the figshare repository

To run the clinical trials pipeline:
    $ kedro run --pipeline data_collection_clinical_trials

To run the iCite data collection pipeline:
    $ kedro run --pipeline data_collection_icite
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes.clinicaltrials import collect_and_save_clinical_trials
from .nodes.icite import collect_and_save_icite_data


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    """
    Create the PubMed and iCite data collection pipeline.

    This pipeline includes both clinical trials data collection and iCite data collection.
    """
    return pipeline(
        [
            # Clinical trials data collection node
            node(
                func=collect_and_save_clinical_trials,
                inputs={
                    "api_download_url": "params:clinical_trials.api_download_url",
                    "zip_save_path": "params:clinical_trials.zip_save_path",
                    "extract_path": "params:clinical_trials.extract_path",
                },
                outputs="clinical_trials_raw",
                name="collect_clinical_trials",
            ),
            # iCite data collection node
            node(
                func=collect_and_save_icite_data,
                inputs={
                    "figshare_download_url": "params:icite.figshare_download_url",
                    "zip_save_path": "params:icite.zip_save_path",
                    "extract_path": "params:icite.extract_path",
                },
                outputs="icite.raw",
                name="collect_icite_data",
            ),
        ],
        namespace="pubmed.data_collection",
    )
