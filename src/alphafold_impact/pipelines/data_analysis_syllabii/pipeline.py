"""
This is a boilerplate pipeline 'data_analysis_syllabii'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    compute_trends,
    plot_trends
)


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([
        node(
            compute_trends,
            inputs="os.labelled_biology_courses.primary",
            outputs=["os.data_analysis.outputs.trends", "os.data_analysis.outputs.regional_trends"],
            tags=["trends"]
        )
    ])
