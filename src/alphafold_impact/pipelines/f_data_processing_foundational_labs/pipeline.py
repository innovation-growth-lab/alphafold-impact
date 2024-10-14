"""
This is a boilerplate pipeline 'data_analysis_lab_productivity'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    get_foundational_lab_staggered_outputs,
)


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:

    staggered_pipeline = pipeline(
        [
            node(
                get_foundational_lab_staggered_outputs,
                inputs={
                    "data_loaders": "foundational_lab.data_collection.publications.raw",
                    "mapping_df": "foundational_lab.data_collection.candidates.map",
                    "publications_data": "publications.data.outputs",
                    "pdb_submissions": "pdb.entries.intermediate",
                    "patents_data": "lens.data_processing.primary",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "institutional_data": "foundational_lab.data_collection.institution_info.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs=[
                    "foundational_lab.data_analysis.staggered.outputs.primary",
                    "foundational_lab.data_analysis.staggered.outputs.quarterly.primary",
                ],
                name="create_foundational_lab_outputs",
            ),
        ],
        tags=["data_processing_foundational_labs"],
    )

    return staggered_pipeline
