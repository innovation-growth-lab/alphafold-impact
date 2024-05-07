"""Pipeline for data analysis.

This pipeline finds Structural and Experimental Biology
NIH funding opportunities.
To run this pipeline, use the following command:

    $ kedro run --pipeline data_analysis_nih_funding_opportunities
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import find_relevant_funding_opps


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=find_relevant_funding_opps,
                inputs=[
                    "nih.data_processing.funding_opportunities.intermediate",
                    "nih.data_processing.funding_opportunities_mesh_tagging_output_df.intermediate",
                    "nih.data_collection.mesh_terms",
                    "params:nih.data_analysis.funding_opportunities.subfield_label_structural_bio",
                    "params:nih.data_analysis.funding_opportunities.mesh_tree_structural_bio",
                    "params:nih.data_analysis.funding_opportunities.search_terms_structural_bio",
                    "params:nih.data_analysis.funding_opportunities.search_cols",
                ],
                outputs="nih.data_analysis.funding_opportunities.structural_biology",
            ),
            node(
                func=find_relevant_funding_opps,
                inputs=[
                    "nih.data_processing.funding_opportunities.intermediate",
                    "nih.data_processing.funding_opportunities_mesh_tagging_output_df.intermediate",
                    "nih.data_collection.mesh_terms",
                    "params:nih.data_analysis.funding_opportunities.subfield_label_experimental_bio",
                    "params:nih.data_analysis.funding_opportunities.mesh_tree_experimental_bio",
                    "params:nih.data_analysis.funding_opportunities.search_terms_experimental_bio",
                    "params:nih.data_analysis.funding_opportunities.search_cols",
                ],
                outputs="nih.data_analysis.funding_opportunities.experimental_biology",
            ),
        ],
    )
