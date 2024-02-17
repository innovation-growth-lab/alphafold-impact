"""Pipeline for data collection.

This pipeline processes NIH funding opportunities data.
To run this pipeline, use the following command:

    $ kedro run --pipeline data_processing_nih_funding_opportunities
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import clean_nih_funding_opportunites

# from ...utils.generic_batch_file import


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=clean_nih_funding_opportunites,
                inputs=[
                    "nih.data_collection.funding_opportunities.raw",
                    "params:nih.data_processing.funding_opportunities.year",
                ],
                outputs="nih.data_processing.funding_opportunities.intermediate",
            ),
        ],
    )
