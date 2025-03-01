"""
This is a boilerplate pipeline 'data_analysis_lab_productivity'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    get_lab_individual_outputs,
)
from ..g_data_collection_authors.nodes import (  # pylint: disable=E0402
    aggregate_to_quarterly,
)


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:

    staggered_pipeline = pipeline(
        [
            node(
                get_lab_individual_outputs,
                inputs={
                    "data_loaders": "foundational_lab.data_collection.publications.raw",
                    "publications_data": "publications.data.outputs",
                    "pdb_submissions": "pdb.entries.primary",
                    "patents_data": "lens.data_processing.primary",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "institutional_data": "foundational_lab.data_collection.institution_info.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs=[
                    "foundational_lab.outputs.primary",
                ],
                name="create_foundational_lab_outputs",
            ),
            node(
                func=aggregate_to_quarterly,
                inputs={
                    "data": "foundational_lab.outputs.primary",
                },
                outputs="foundational_lab.outputs.quarterly",
                name="aggregate_foundational_lab_to_quarterly",
            ),
        ],
        tags=["data_processing_foundational_labs"],
    )

    return staggered_pipeline
