"""Pipeline for data collection.

This pipeline processes NIH funding opportunities data.
To run this pipeline, use the following command:

    $ kedro run --pipeline data_processing_nih_funding_opportunities
"""

from kedro.pipeline import Pipeline, pipeline, node
from alphafold_impact.utils.nih_mesh_tagging import (
    generate_text_for_mesh_tagging,
    skr_web_python_api_generic_batch,
    mesh_results_to_df,
)
from .nodes import clean_nih_funding_opportunites


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
            node(
                func=generate_text_for_mesh_tagging,
                inputs=[
                    "nih.data_processing.funding_opportunities.intermediate",
                    "params:nih.data_processing.funding_opportunities.mesh_tagging.id_col",
                    "params:nih.data_processing.funding_opportunities.mesh_tagging.text_col",
                ],
                outputs="nih.data_processing.funding_opportunities_mesh_tagging_input_text.intermediate",
            ),
            node(
                func=skr_web_python_api_generic_batch,
                inputs=[
                    "params:mesh_tagging.email",
                    "params:mesh_tagging.umls_api_key",
                    "params:nih.data_processing.funding_opportunities.mesh_tagging.batch_name",
                    "nih.data_processing.funding_opportunities_mesh_tagging_input_text.intermediate",
                    "params:nih.data_processing.funding_opportunities.mesh_tagging.cmd",
                    "params:nih.data_processing.funding_opportunities.mesh_tagging.cmdargs",
                ],
                outputs="nih.data_processing.funding_opportunities_mesh_tagging_output_text.intermediate",
            ),
            node(
                func=mesh_results_to_df,
                inputs=[
                    "nih.data_processing.funding_opportunities_mesh_tagging_output_text.intermediate",
                ],
                outputs="nih.data_processing.funding_opportunities_mesh_tagging_output_df.intermediate",
            ),
        ],
    )
