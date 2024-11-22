"""
This is a boilerplate pipeline 'data_processing_pdb'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import collect_pdb_details, merge_uniprot_data


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    return pipeline(
        [
            node(
                collect_pdb_details,
                inputs={
                    "pdb_df": "pdb.entries.raw",
                    "api_config": "params:oa.data_collection.depth.api",
                },
                outputs="pdb.entries.intermediate",
                name="fetch_oa_details_from_pdb_pubs",
            ),
            node(
                merge_uniprot_data,
                inputs={
                    "pdb_df": "pdb.entries.intermediate",
                    "uniprot_df": "uniprot.entries.raw",
                },
                outputs="pdb.entries.processed",
                name="merge_uniprot_data",
            ),
        ],
        tags=["data_processing_pdb"],
    )
