
from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    get_applied_lab_staggered_outputs,
)


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:

    staggered_pipeline = pipeline(
        [
            node(
                get_applied_lab_staggered_outputs,
                inputs={
                    "data_loaders": "applied_lab.data_collection.publications.raw",
                    "mapping_df": "applied_lab.data_collection.candidates.map",
                    "publications_data": "publications.data.outputs",
                    "pdb_submissions": "pdb.entries.primary",
                    "patents_data": "lens.data_processing.primary",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "institutional_data": "applied_lab.data_collection.institution_info.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs=[
                    "applied_lab.data_analysis.staggered.outputs.primary",
                    "applied_lab.data_analysis.staggered.outputs.quarterly.primary",
                ],
                name="create_applied_lab_outputs",
            ),
        ],
        tags=["data_processing_applied_labs"],
    )

    return staggered_pipeline