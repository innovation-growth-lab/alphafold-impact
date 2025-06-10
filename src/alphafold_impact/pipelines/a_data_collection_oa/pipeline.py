"""
Pipeline for collecting OpenAlex data. It includes:

    Base:
    - Collecting OpenAlex publications given a work_id.
    - Collect references and citations for a work_id.
    - Collect citations for all citations of a work_id.

    GtR:
        - Collecting OpenAlex publications that match the doi values for the Gateway to Research
            publications.
        - Collecting all papers that these publications cite.

PIPELINE FLOW:
    STEP 1: This pipeline (a_data_collection_oa) collects raw OpenAlex data
    NEXT: → b__data_processing_oa (processes the raw OpenAlex data)
    THEN: → b_data_processing_baselines (creates counterfactual candidates)
    THEN: → BACK to a_data_collection_oa (collects more data for counterfactuals)
    THEN: → b__data_processing_oa (processes new counterfactual data)
    THEN: → c_data_collection_s2 (adds citation intent from Semantic Scholar)
    FINALLY: → d_data_processing_chains (final citation chain analysis)

To run this pipeline, use the following command:

    $ kedro run --pipeline data_collection_oa

To run a specific pipeline, use a combination of the --namespace and --tags flags:

    $ kedro run --pipeline data_collection_oa --namespace oa.data_collection.gtr
    $ kedro run --pipeline data_collection_oa --tags downstream_impact



"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import (
    fetch_subfield_baseline,
    fetch_citation_to_specific_depth,
    preprocess_baseline_data,
    preprocess_ct_level1_data,
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613

    # STEP 1A: Collect AlphaFold citations to specific depth
    # → NEXT: Goes to b__data_processing_oa for processing
    fixed_depth_alphafold_pipeline = pipeline(
        [
            node(
                func=fetch_citation_to_specific_depth,
                inputs={
                    "seed_paper": "params:oa.data_collection.depth.get.work_ids",
                    "api_config": "params:oa.data_collection.depth.api",
                    "filter_config": "params:oa.data_collection.depth.filter",
                    "max_depth": "params:oa.data_collection.depth.af.max_depth",
                },
                outputs="oa.data_collection.triad.depth.level.raw",
                name="fetch_af_citation_to_specific_depth",
            )
        ],
        tags="data_collection_oa",
    )

    # STEP 1B: Collect structural biology baseline papers
    # → NEXT: Goes to b__data_processing_oa for processing
    structural_biology_pipeline = pipeline(
        [
            node(
                func=fetch_subfield_baseline,
                inputs={
                    "oa_concept_ids": "params:oa.data_collection.subfield.structural_biology.concept_ids",  # pylint: disable=line-too-long
                    "from_publication_date": "params:oa.data_collection.subfield.structural_biology.from_date",  # pylint: disable=line-too-long
                    "api_config": "params:oa.data_collection.subfield.structural_biology.api",
                },
                outputs="oa.data_collection.subfield.structural_biology.raw",
                name="fetch_baseline_sb",
            )
        ],
        tags="data_collection_oa",
    )

    # STEP 1C: Collect citations for structural biology baseline papers
    # → NEXT: Goes to b__data_processing_oa for processing
    fixed_depth_structural_biology_pipeline = pipeline(
        [
            node(
                func=preprocess_baseline_data,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.raw",
                    "alphafold_papers": "params:oa.data_collection.depth.get.work_ids",
                },
                outputs="oa.structural_biology.seed_papers",
                name="preprocess_sb",
            ),
            node(
                func=fetch_citation_to_specific_depth,
                inputs={
                    "seed_paper": "oa.structural_biology.seed_papers",
                    "api_config": "params:oa.data_collection.depth.api",
                    "filter_config": "params:oa.data_collection.depth.filter",
                    "max_depth": "params:oa.data_collection.depth.sb.max_depth",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.raw",
                name="fetch_sb_citation_to_specific_depth",
            ),
        ],
        tags=["data_collection_oa"],
    )

    # STEP 4: Collect additional citations for CT papers (runs AFTER b_data_processing_baselines)
    # This step is triggered by data created in b_data_processing_baselines
    # → NEXT: Goes back to b__data_processing_oa for processing the new counterfactual data
    last_level_ct_pipeline = pipeline(
        [
            node(
                func=preprocess_ct_level1_data,
                inputs={
                    "data": "oa.data_processing.structural_biology.depth.reassigned.ct.intermediate",
                    "alphafold_papers": "params:oa.data_collection.depth.get.work_ids",
                },
                outputs=["oa.ct.level1.seed_papers", "oa.ct.level1.seen_papers"],
                name="preprocess_ct_level1",
            ),
            node(
                func=fetch_citation_to_specific_depth,
                inputs={
                    "seed_paper": "oa.ct.level1.seed_papers",
                    "papers_seen": "oa.ct.level1.seen_papers",
                    "api_config": "params:oa.data_collection.depth.api",
                    "filter_config": "params:oa.data_collection.depth.filter",
                    "start_level": "params:oa.data_collection.depth.levels.2",
                    "max_depth": "params:oa.data_collection.depth.af.max_depth",
                },
                outputs="oa.data_collection.structural_biology.depth.ct.2.ptd.raw",
                name="fetch_ct_level1_citation_to_specific_depth",
            ),
        ],
        tags=["last_level_ct_pipeline"],
    )

    return (
        fixed_depth_alphafold_pipeline
        + structural_biology_pipeline  # SB pipeline
        + fixed_depth_structural_biology_pipeline
        + last_level_ct_pipeline
    )
