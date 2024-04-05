"""
This is a boilerplate pipeline 'data_collection_s2'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import (
    get_citation_intent_from_oa_dataset
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613

    primary_data_intent_pipeline = pipeline(
        [
            node(
                func=get_citation_intent_from_oa_dataset,
                inputs={
                    "oa_dataset": "oa.data_processing.depth.intermediate",
                    "base_url": "params:s2.data_collection.strength.api.base_url",
                    "fields": "params:s2.data_collection.strength.api.fields",
                    "api_key": "params:s2.data_collection.strength.api.key",
                    "perpage": "params:s2.data_collection.strength.api.perpage",
                },
                outputs="oa.data_processing.depth.primary",
            ),
        ],
        tags=["s2.primary.intent"],
    )



    return primary_data_intent_pipeline
