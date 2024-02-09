"""Pipeline for data collection.

This pipeline fetches NIH funding opportunities data.
To run this pipeline, use the following command:

    $ kedro run --pipeline data_collection_nih_funding_opportunities
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import collect_nih_funding_opportunities


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=collect_nih_funding_opportunities,
                inputs={"api_url": "params:api_url", "end_date": "params:end_date"},
                outputs="raw",
            ),
        ],
        namespace="nih.data_collection.funding_opportunities",
    )
