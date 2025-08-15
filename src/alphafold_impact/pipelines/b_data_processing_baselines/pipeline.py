"""
The poorly named pipeline for processing the structural biology baseline data to extract
the counterfactual papers to AlphaFold.

PIPELINE FLOW:
    STEP 4: This pipeline (b_data_processing_baselines) creates counterfactual candidates
    FROM: ← b__data_processing_oa (provides processed OpenAlex data)
    NEXT: → a_data_collection_oa (triggers collection of more citation data for counterfactuals)
    THEN: → b__data_processing_oa (processes the new counterfactual data)
    THEN: → c_data_collection_s2 (adds citation intent for counterfactuals)
    FINALLY: → d_data_processing_chains (final citation chain analysis)
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import process_af_data, process_baseline_data, assign_focal_label


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    return pipeline(
        [
            # STEP 4A: Process baseline structural biology data to find counterfactual candidates
            # ← FROM: b__data_processing_oa (provides processed baseline and AlphaFold data)
            node(
                process_baseline_data,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.af.primary",
                    "intent_data": "oa.data_processing.depth.level.0.baseline.primary",
                    "baseline_data": "oa.data_collection.subfield.structural_biology.depth.0.intermediate", # actually level 1 (eventually)
                    "seed_baseline_data": "oa.data_processing.subfield.structural_biology.primary",
                },
                outputs="baseline_candidates",
                name="process_sb_for_chain_assignment",
            ),
            # STEP 4B: Process AlphaFold data to create target profile for comparison
            # ← FROM: b__data_processing_oa (provides processed AlphaFold data)
            node(
                process_af_data,
                inputs={"alphafold_data": "oa.data_processing.depth.af.primary"},
                outputs="alphafold_target",
                name="process_af_for_chain_assignment",
            ),
            # STEP 4C: Assign labels and select best counterfactual papers
            # → NEXT: Triggers a_data_collection_oa to collect more citation data for these papers
            node(
                assign_focal_label,
                inputs={
                    "baseline_candidates": "baseline_candidates",
                    "alphafold_target": "alphafold_target",
                },
                outputs="chains.seed_technologies.intermediate",
                name="assign_ct_focal_label",
            ),
        ],
        tags=["data_processing_baselines"],
    )
