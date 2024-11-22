"""
This is a boilerplate pipeline 'data_collection_pdb'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import fetch_uniprot_details


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    return pipeline(
        [
            node(
                func=fetch_uniprot_details,
                inputs={
                    "pdb_data": "pdb.entries.raw",
                    },
                outputs="uniprot.entries.raw",
                name="fetch_uniprot_details",
            )
        ],
        tags=["data_collection_uniprot"],
    )
