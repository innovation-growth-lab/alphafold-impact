"""Pipeline for ClinicalTrials.gov data collection.

This pipeline fetches ClinicalTrials.gov data using the new API v2 endpoints.
The legacy endpoints are also supported for backward compatibility.

To run this pipeline, use the following command:
    $ kedro run --pipeline data_collection_clinical_trials

Migration Notes:
- The new API v2 endpoint provides JSON format data with a normalized schema
- Legacy XML endpoints are still supported but deprecated
- Recommended endpoint: /api/v2/studies/download?format=json.zip
- Legacy endpoint: /AllAPIXML.zip (still works but deprecated)
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    collect_and_save_clinical_trials_v2,
    collect_and_save_nih_clinical_trials  # Keep for backward compatibility
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    """
    Create the ClinicalTrials.gov data collection pipeline.
    
    This pipeline supports both the new API v2 and legacy endpoints.
    The function used depends on the parameters provided.
    """
    return pipeline(
        [
            node(
                func=collect_and_save_clinical_trials_v2,
                inputs={
                    "api_v2_download_url": "params:api_v2_download_url",
                    "zip_save_path": "params:zip_save_path", 
                    "extract_path": "params:extract_path",
                },
                outputs="raw_v2",
                name="collect_clinical_trials_v2",
            ),
            # Legacy node for backward compatibility
            node(
                func=collect_and_save_nih_clinical_trials,
                inputs={
                    "nih_ct_download_url": "params:download_url",
                    "nih_ct_zip_save_path": "params:zip_save_path",
                    "nih_ct_extract_path": "params:extract_path",
                },
                outputs="raw_legacy",
                name="collect_clinical_trials_legacy",
            ),
        ],
        namespace="pubmed.data_collection.icite",
    )
