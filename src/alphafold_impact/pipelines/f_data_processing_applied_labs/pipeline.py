from kedro.pipeline import Pipeline, pipeline, node
from ..f_data_processing_foundational_labs.nodes import (  # pylint: disable=E0402
    get_lab_individual_outputs,
    aggregate_to_quarterly
)


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:

    staggered_pipeline = pipeline(
        [
            node(
                get_lab_individual_outputs,
                inputs={
                    "data_loaders": "applied_lab.data_collection.publications.raw",
                    "publications_data": "publications.data.outputs",
                    "pdb_submissions": "pdb.entries.primary",
                    "patents_data": "lens.data_processing.primary",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "institutional_data": "applied_lab.data_collection.institution_info.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="applied_lab.outputs.primary",
                name="create_applied_lab_outputs",
            ),
            node(
                func=aggregate_to_quarterly,
                inputs={
                    "data": "applied_lab.outputs.primary",
                },
                outputs="applied_lab.outputs.quarterly",
                name="aggregate_applied_lab_to_quarterly",
            ),
        ],
        tags=["data_processing_applied_labs"],
    )

    return staggered_pipeline
