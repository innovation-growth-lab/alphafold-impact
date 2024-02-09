"""Pipeline for data collection.

This pipeline processes NIH funding opportunities data.
To run this pipeline, use the following command:

    $ kedro run --pipeline data_processing_nih_funding_opportunities
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import clean_nih_funding_opportunites


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=clean_nih_funding_opportunites,
                inputs=None,
                outputs="raw",
            ),
        ],
        namespace="nih.data_processing.funding_opportunities",
    )
