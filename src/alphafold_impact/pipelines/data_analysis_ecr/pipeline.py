"""
This is a boilerplate pipeline 'data_analysis_ecr'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import run_analysis


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                run_analysis,
                inputs={
                    "data_loaders": "ecr.authors.publications.raw",
                    "af_data": "analysis.af.level0_data",
                },
                outputs=[
                    "ecr.data_analysis.subfield.af",
                    "ecr.data_analysis.subfield.no_af",
                    "ecr.data_analysis.field.af",
                    "ecr.data_analysis.field.no_af",
                    "ecr.data_analysis.counts.af",
                    "ecr.data_analysis.counts.no_af",
                    "ecr.data_analysis.citations.af",
                    "ecr.data_analysis.citations.no_af",
                ],
                tags=["ecr_concatenate_data", "ecr_analysis"],
            ),
        ]
    )
