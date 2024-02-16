"""Pipeline for data collection.

This pipeline fetches NIH funded projects data.
To run this pipeline, use the following command:

    $ kedro run --pipeline data_collection_nih_funded_projects
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import collect_nih_funded_projects


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=collect_nih_funded_projects,
                inputs={"years": "params:years", "api_url": "params:api_url"},
                outputs="raw",
            ),
        ],
        namespace="nih.data_collection.funded_projects",
    )
