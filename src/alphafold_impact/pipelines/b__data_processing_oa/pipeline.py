"""This module contains the pipeline for data processing for OA data.

The pipeline contains nodes for processing data by level, combining data
from different levels, and mesh tagging.

PIPELINE FLOW:
    STEP 2: This pipeline (b__data_processing_oa) processes raw OpenAlex data
    FROM: ← a_data_collection_oa (provides raw OpenAlex data)
    NEXT: → b_data_processing_baselines (creates counterfactual candidates)
    ALSO NEXT: → c_data_collection_s2 (adds citation intent from Semantic Scholar)
    LATER: ← a_data_collection_oa (provides additional counterfactual data)
    THEN: → c_data_collection_s2 (processes counterfactual citation intent)
    FINALLY: → d_data_processing_chains (final citation chain analysis)

To run this pipeline, use the following command:

    $ kedro run --pipeline oa.data_processing.depth
"""

from kedro.pipeline import Pipeline, node, pipeline
from alphafold_impact import settings  # pylint: disable=import-error

from .nodes import (
    process_subfield_data,
    process_data_by_level,
    combine_levels_data,
    combine_levels_data_counterfactuals,
    reassign_ct_levels,
    process_data_by_level_ptd,
    concat_pq_ptd,
)


def create_pipeline(  # pylint: disable=unused-argument&missing-function-docstring
    **kwargs,
) -> Pipeline:

    # STEP 2A: Process AlphaFold citation data by level
    # ← FROM: a_data_collection_oa (oa.data_collection.triad.depth.level.raw)
    # → NEXT: Combined in 2B and sent to c_data_collection_s2 for citation intent analysis
    af_processing_pipeline = pipeline(
        [
            node(
                func=process_data_by_level,
                inputs={
                    "data": "oa.data_collection.triad.depth.level.raw",
                    "level": f"params:oa.data_collection.depth.levels.{level}",
                },
                outputs=f"oa.data_processing.depth.{level}.intermediate",
                name=f"process_af_level_{str(level)}",
            )
            for level in settings.DYNAMIC_PIPELINES_MAPPING["depth_levels"]
        ],
        tags=["data_processing_oa", "Q1"],
    )

    # STEP 2B: Combine AlphaFold levels for processing
    # → NEXT: Goes to c_data_collection_s2 for citation intent analysis
    combine_levels_pipeline = pipeline(
        [
            node(
                func=combine_levels_data,
                inputs={
                    # "unique": "params:false_",
                    "level0": "oa.data_processing.depth.0.intermediate",
                    "level1": "oa.data_processing.depth.1.intermediate",
                    "level2": "oa.data_processing.depth.2.intermediate",
                },
                outputs="oa.data_processing.depth.intermediate",
                name="combine_af_levels",
            )
        ],
        tags=["data_processing_oa", "Q1"],
    )

    # STEP 2C: Process structural biology baseline data
    # ← FROM: a_data_collection_oa (oa.data_collection.subfield.structural_biology.raw)
    # → NEXT: Goes to b_data_processing_baselines and c_data_collection_s2
    structural_biology_processing_pieline = pipeline(
        [
            node(
                func=process_subfield_data,
                inputs=["oa.data_collection.subfield.structural_biology.raw"],
                outputs="oa.data_processing.subfield.structural_biology.primary",
                name="process_sb",
            )
        ],
        tags=["data_processing_oa"],
    )

    # STEP 2D: Process baseline level 0 data
    # ← FROM: a_data_collection_oa
    # → NEXT: Combined with other levels and sent to b_data_processing_baselines
    baseline_level0_pipeline = pipeline(
        [
            node(
                func=process_data_by_level,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.raw",
                    "level": "params:oa.data_collection.depth.levels.0",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.0.intermediate",
                name="process_sb_level_0",
            )
        ],
        tags=["data_processing_oa"],
    )

    # STEP 2E: Process baseline level 1 data (large dataset, requires partitioning)
    # ← FROM: a_data_collection_oa
    # → NEXT: Combined with other levels and sent to b_data_processing_baselines
    baseline_level1_pipeline = pipeline(
        [
            node(
                func=process_data_by_level_ptd,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.raw",
                    "level": "params:oa.data_collection.depth.levels.1",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.1.ptd.intermediate",
                name="process_sb_level_1",
            ),
            node(
                func=concat_pq_ptd,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.depth.1.ptd.intermediate",  # pylint: disable=line-too-long
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.1.intermediate",
                name="concatenate_sb_partitioned",
            ),
        ],
        tags=[
            "data_processing_oa",
        ],
    )

    # STEP 2F: Combine baseline levels and split into counterfactual vs other papers
    # → NEXT: Goes to b_data_processing_baselines (for counterfactual identification)
    # → ALSO: The output is used by a_data_collection_oa for additional data collection
    reassign_baseline_levels_pipeline = pipeline(
        [
            node(
                func=combine_levels_data,
                inputs={
                    "level_seed": "oa.data_processing.subfield.structural_biology.primary",
                    "level0": "oa.data_collection.subfield.structural_biology.depth.0.intermediate",
                    "level1": "oa.data_collection.subfield.structural_biology.depth.1.intermediate",
                },
                outputs="oa.data_processing.structural_biology.depth.intermediate",
                name="combine_sb_levels",
            ),
            node(
                func=reassign_ct_levels,
                inputs={
                    "data": "oa.data_processing.structural_biology.depth.intermediate",
                    "ct_data": "chains.seed_technologies.intermediate",
                },
                outputs=[
                    "oa.data_processing.structural_biology.depth.reassigned.ct.intermediate",
                    "oa.data_processing.structural_biology.depth.other.intermediate",
                ],
                name="extract_ct_sb_levels",
            ),
        ],
        tags=["reassign_ct_levels", "ct_collection"],
    )

    # STEP 5: Process additional CT data (runs AFTER a_data_collection_oa collects more data)
    # ← FROM: a_data_collection_oa (oa.data_collection.structural_biology.depth.ct.2.ptd.raw)
    # → NEXT: Goes to c_data_collection_s2 for citation intent analysis
    post_level2_dw_pipeline = pipeline(
        [
            node(
                func=process_data_by_level_ptd,
                inputs={
                    "data": "oa.data_collection.structural_biology.depth.ct.2.ptd.raw",
                    "level": "params:oa.data_collection.depth.levels.2",
                    # HACK: Legacy from when SB had papers to lvl 3, current input only has lvl 2
                },
                outputs="oa.data_collection.structural_biology.depth.ct.2.ptd.intermediate",
                name="process_ct_sb_level_2",
            ),
            node(
                func=concat_pq_ptd,
                inputs={
                    "data": "oa.data_collection.structural_biology.depth.ct.2.ptd.intermediate",
                },
                outputs="oa.data_collection.structural_biology.depth.ct.2.intermediate",
                name="concatenate_ct_sb_partitioned",
            ),
            node(
                func=combine_levels_data_counterfactuals,
                inputs={
                    "level1": "oa.data_processing.structural_biology.depth.reassigned.ct.intermediate",  # pylint: disable=line-too-long
                    "level2": "oa.data_collection.structural_biology.depth.ct.2.intermediate",
                },
                outputs="oa.data_processing.structural_biology.depth.ct.intermediate",
                name="combine_ct_sb_levels",
            ),
        ],
        tags=["collect_new_level_for_ct", "ct_collection"],
    )

    return (
        af_processing_pipeline
        + combine_levels_pipeline
        + structural_biology_processing_pieline
        + baseline_level0_pipeline
        + baseline_level1_pipeline
        + reassign_baseline_levels_pipeline
        + post_level2_dw_pipeline
    )
