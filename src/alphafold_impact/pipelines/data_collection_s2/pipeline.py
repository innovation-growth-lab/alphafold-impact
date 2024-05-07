"""This module contains the Kedro pipeline for the data collection from Semantic Scholar API.

The pipeline contains a single node that retrieves the citation intent from the OpenAlex
dataset.

To run this pipeline, use the following command:

    $ kedro run --pipeline s2.primary.intent

"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import (
    get_citation_intent_from_oa_dataset,
    get_baseline_seed_intent
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
                outputs="oa.data_processing.depth.all.primary",
            ),
        ],
        tags=["s2.intent.primary"],
    )

    baseline_seed_pipeline = pipeline(
        [
            node(
                func=get_baseline_seed_intent,
                inputs={
                    "oa_dataset": "oa.data_processing.subfield.structural_biology.primary",
                    "base_url": "params:s2.data_collection.strength.api.base_url",
                    "fields": "params:s2.data_collection.strength.api.fields",
                    "api_key": "params:s2.data_collection.strength.api.key",
                    "perpage": "params:s2.data_collection.strength.api.perpage",
                },
                outputs="oa.data_processing.depth.level.0.baseline.primary",
            ),
        ],
        tags=["s2.intent.baseline"],
    )

    other_sb_data_intent_pipeline = pipeline(
        [
            node(
                func=get_citation_intent_from_oa_dataset,
                inputs={
                    "oa_dataset": "oa.data_processing.structural_biology.depth.reassigned.other.intermediate",
                    "base_url": "params:s2.data_collection.strength.api.base_url",
                    "fields": "params:s2.data_collection.strength.api.fields",
                    "api_key": "params:s2.data_collection.strength.api.key",
                    "perpage": "params:s2.data_collection.strength.api.perpage",
                },
                outputs="oa.data_processing.depth.other.primary",
            ),
        ],
        tags=["s2_intent_other_sb"],
    )

    return primary_data_intent_pipeline + baseline_seed_pipeline + other_sb_data_intent_pipeline
