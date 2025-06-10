"""This module contains the Kedro pipeline for the data collection from Semantic Scholar API.

The pipeline contains nodes that retrieve citation intent from the OpenAlex dataset using
Semantic Scholar's API to understand how papers cite each other (background, methodology, result).

PIPELINE FLOW:
    STEP 6: This pipeline (c_data_collection_s2) adds citation intent data from Semantic Scholar
    FROM: ← b__data_processing_oa (provides processed OpenAlex data)
    NEXT: → d_data_processing_chains (final citation chain analysis)

    The pipeline runs multiple times:
    - First for AlphaFold papers (after step 2)
    - Then for baseline seed papers (after step 2)
    - Finally for counterfactual and other papers (after step 5)

To run this pipeline, use the following command:

    $ kedro run --pipeline s2.primary.intent

"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import (
    get_citation_intent_from_oa_dataset,
    process_citation_levels,
    get_baseline_seed_intent,
    get_baseline_intent,
)
from .nodes_bulk import (
    get_s2_presigned_urls,
    get_s2_bulk_data,
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613

    # STEP 6A: Add citation intent data for AlphaFold papers
    # ← FROM: b__data_processing_oa (oa.data_processing.depth.intermediate)
    # → NEXT: Goes to d_data_processing_chains for final analysis
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
                outputs="s2.af.intents.intermediate",
                name="fetch_s2_citation_intent_af",
            ),
            node(
                func=process_citation_levels,
                inputs={
                    "oa_dataset": "oa.data_processing.depth.intermediate",
                    "levels": "s2.af.intents.intermediate",
                },
                outputs="oa.data_processing.depth.af.primary",
                name="process_s2_citation_intent_af",
            ),
        ],
        tags=["data_collection_s2"],
    )

    # STEP 6B: Add citation intent data for baseline seed papers
    # ← FROM: b__data_processing_oa (oa.data_processing.subfield.structural_biology.primary)
    # → NEXT: Used by b_data_processing_baselines for counterfactual selection
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
                name="fetch_s2_citation_intent_baseline_seed",
            )
        ],
        tags=["data_collection_s2_baseline_seed"],
    )

    # STEP 6C: Add citation intent data for counterfactual papers
    # ← FROM: b__data_processing_oa (oa.data_processing.structural_biology.depth.ct.intermediate)
    # → NEXT: Goes to d_data_processing_chains for final analysis
    counterfactual_pipeline = pipeline(
        [
            node(
                func=get_baseline_intent,
                inputs={
                    "oa_dataset": "oa.data_processing.structural_biology.depth.ct.intermediate",
                    "base_url": "params:s2.data_collection.strength.api.base_url",
                    "fields": "params:s2.data_collection.strength.api.fields",
                    "api_key": "params:s2.data_collection.strength.api.key",
                    "perpage": "params:s2.data_collection.strength.api.perpage",
                },
                outputs="s2.baseline.ct.intents.intermediate",
                name="fetch_s2_citation_intent_baseline_ct",
            ),
            node(
                func=process_citation_levels,
                inputs={
                    "oa_dataset": "oa.data_processing.structural_biology.depth.ct.intermediate",
                    "levels": "s2.baseline.ct.intents.intermediate",
                },
                outputs="oa.data_processing.depth.ct.primary",
                name="process_s2_citation_intent_baseline_ct",
            ),
        ],
        tags=["ct_s2_pipeline"],
    )

    # STEP 6D: Add citation intent data for other structural biology papers
    # ← FROM: b__data_processing_oa (oa.data_processing.structural_biology.depth.other.intermediate)
    # → NEXT: Goes to d_data_processing_chains for final analysis
    other_pipeline = pipeline(
        [
            node(
                func=get_baseline_intent,
                inputs={
                    "oa_dataset": "oa.data_processing.structural_biology.depth.other.intermediate",
                    "base_url": "params:s2.data_collection.strength.api.base_url",
                    "fields": "params:s2.data_collection.strength.api.fields",
                    "api_key": "params:s2.data_collection.strength.api.key",
                    "perpage": "params:s2.data_collection.strength.api.perpage",
                },
                outputs="s2.baseline.other.intents.intermediate",
                name="fetch_s2_citation_intent_baseline_other",
            ),
            node(
                func=process_citation_levels,
                inputs={
                    "oa_dataset": "oa.data_processing.structural_biology.depth.other.intermediate",
                    "levels": "s2.baseline.other.intents.intermediate",
                },
                outputs="oa.data_processing.depth.other.primary",
                name="process_s2_citation_intent_baseline_other",
            ),
        ],
        tags=["sb_s2_pipeline"],
    )

    # STEP 6E: Alternative bulk data collection method (optional)
    # This provides an alternative way to collect Semantic Scholar data in bulk
    bulk_data_pipeline = pipeline(
        [
            node(
                func=get_s2_presigned_urls,
                inputs={
                    "bulk_url": "params:s2.data_collection.strength.api.bulk_url",
                    "api_key": "params:s2.data_collection.strength.api.key",
                },
                outputs="s2.data_collection.strength.api.bulk_urls",
            ),
            node(
                func=get_s2_bulk_data,
                inputs="s2.data_collection.strength.api.bulk_urls",
                outputs="s2.data_collection.strength.api.bulk_data",
                name="fetch_s2_bulk_data",
            ),
        ],
        tags=["bulk_data_s2"],
    )
    return (
        primary_data_intent_pipeline
        + baseline_seed_pipeline
        + counterfactual_pipeline
        + other_pipeline
        + bulk_data_pipeline
    )
