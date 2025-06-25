"""
This is a boilerplate pipeline 'publications_descriptive_translational'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    generate_fig_counts_and_field_shares,
    generate_fig_all_counts,
    generate_fig_researcher_counts,
    generate_fig_within_plots,
    generate_fig3,
    generate_summary_tables_latex,
    plot_regression_results,
    combined_publications_researchers_pp,
    generate_fig_descriptive_protein_charspy,
    generate_fig_descriptive_translational,
    generate_fig_topics_over_time
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    create_publications_data_pipeline = pipeline(
        [
            node(
                func=generate_fig_counts_and_field_shares, # current Fig 4
                inputs={"publications": "publications.data.outputs"},
                outputs="fig.counts_and_field_shares",
                name="generate_fig_counts_and_field_shares",
            ),
            node(
                func=generate_fig_all_counts,
                inputs={"publications": "publications.data.outputs"},
                outputs="fig.all_counts",
                name="generate_fig_all_counts",
            ),
            node(
                func=generate_fig_researcher_counts,
                inputs={"publications": "publications.data.outputs"},
                outputs="fig.researcher_counts",
                name="generate_fig_researcher_counts",
            ),
            node(
                func=generate_fig_within_plots,
                inputs={"publications": "publications.data.outputs"},
                outputs="fig.within_plots",
                name="generate_fig_within_plots",
            ),
            node(
                func=combined_publications_researchers_pp,
                inputs={
                    "publications": "publications.data.outputs",
                    "researchers": "nonecr.publications.quarterly",
                },
                outputs="fig.descriptive_pp",
                name="combined_descriptive_pp",
            ),
            node(
                func=generate_fig_descriptive_protein_charspy,
                inputs={"publications": "publications.data.outputs"},
                outputs="fig.descriptive_protein_charspy",
                name="generate_fig_descriptive_protein_charspy",
            ),
            node(
                func=generate_fig_descriptive_translational,
                inputs={"publications": "publications.data.outputs"},
                outputs="fig.descriptive_translational",
                name="generate_fig_descriptive_translational",
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
                    "foundational_labs": "foundational_lab.outputs.quarterly",
                    "applied_labs": "applied_lab.outputs.quarterly",
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
                func=generate_fig_topics_over_time,
                inputs={"publications": "publications.data.outputs"},
                outputs="fig.topics_over_time",
                name="generate_fig_topics_over_time",
            ),
        ]
    )

    return create_publications_data_pipeline
