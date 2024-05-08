"""Pipeline for reporting.

This pipeline plots and saves Funding Opportunity related charts.
To run this pipeline, use the following command:

    $ kedro run --pipeline reporting_nih_funding_opportunities
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    plot_eb_sb_funding_opportunities,
    plot_eb_sb_funding_opportunity_amounts,
    plot_total_funding_opportunities,
    plot_total_funding_opportunity_amounts,
)


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=plot_eb_sb_funding_opportunities,
                inputs=[
                    "nih.data_analysis.funding_opportunities.structural_biology",
                    "nih.data_analysis.funding_opportunities.experimental_biology",
                ],
                outputs="nih.reporting.funding_opportunities.sb_eb",
            ),
            node(
                func=plot_eb_sb_funding_opportunity_amounts,
                inputs=[
                    "nih.data_analysis.funding_opportunities.structural_biology",
                    "nih.data_analysis.funding_opportunities.experimental_biology",
                ],
                outputs="nih.reporting.funding_opportunity_amounts.sb_eb",
            ),
            node(
                func=plot_total_funding_opportunities,
                inputs="nih.data_processing.funding_opportunities.intermediate",
                outputs="nih.reporting.funding_opportunities.total",
            ),
            node(
                func=plot_total_funding_opportunity_amounts,
                inputs="nih.data_processing.funding_opportunities.intermediate",
                outputs="nih.reporting.funding_opportunity_amounts.total",
            ),
        ],
    )
