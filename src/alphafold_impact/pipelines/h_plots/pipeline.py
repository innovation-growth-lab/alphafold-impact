"""
This is a boilerplate pipeline 'publications_descriptive_translational'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    generate_fig1,
    generate_fig2,
    generate_fig3,
    generate_summary_tables_latex,
    plot_regression_results,
    generate_researcher_charts,
    generate_publication_charts,
    combined_publications_researchers_charts
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    create_publications_data_pipeline = pipeline(
        [
            node(
                func=generate_fig1,
                inputs={"publications": "publications.data.outputs"},
                outputs="fig.1",
                name="generate_fig1",
            ),
            node(
                func=generate_fig2,
                inputs={"publications": "publications.data.outputs"},
                outputs="fig.2",
                name="generate_fig2",
            ),
            node(
                func=generate_fig3,
                inputs={
                    "ecrs": "ecr.publications.quarterly",
                    "nonecrs": "nonecr.publications.quarterly",
                },
                outputs="fig.3",
                name="generate_fig3",
            ),
            node(
                func=generate_summary_tables_latex,
                inputs={
                    "publications": "publications.data.outputs",
                    "early_career_researchers": "ecr.publications.quarterly",
                    "established_researchers": "nonecr.publications.quarterly",
                    "foundational_labs": "foundational_lab.data_analysis.staggered.outputs.quarterly.primary",
                    "applied_labs": "applied_lab.data_analysis.staggered.outputs.quarterly.primary",
                },
                outputs="summary_table",
                name="generate_summary_table",
            ),
            node(
                func=plot_regression_results,
                inputs={
                    "reg": "organism_rarity_mean_coef_table.rds",
                    "depth": "params:depth_allgroups",
                    "field": "params:field_all",
                    "subgroup": "params:cem",
                },
                outputs="regression_results",
                name="plot_regression_results",
            ),
            node(
                func=generate_researcher_charts,
                inputs={
                    "ecrs": "ecr.publications.quarterly",
                    "researchers": "nonecr.publications.quarterly"
                },
                outputs="fig.4",
                name="generate_researcher_charts",
            ),
            node(
                func=generate_publication_charts,
                inputs={
                    "publications": "publications.data.outputs",
                },
                outputs="fig.5",
                name="generate_publication_charts",
            ),
            node(
                func=combined_publications_researchers_charts,
                inputs={
                    "publications": "publications.data.outputs",
                    "ecrs": "ecr.publications.quarterly",
                    "researchers": "nonecr.publications.quarterly",
                    "columns": "params:columns.fig_4",
                    "labels": "params:labels.fig_4",
                    "chart_titles": "params:chart_titles.fig_4",
                },
                outputs="fig.4combined",
                name="combined_publications_researchers_charts",
            )
        ]
    )

    return create_publications_data_pipeline
