"""
This is a boilerplate pipeline 'publications_descriptive_translational'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    generate_fig1, generate_fig2
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
        ]
    )

    return create_publications_data_pipeline
