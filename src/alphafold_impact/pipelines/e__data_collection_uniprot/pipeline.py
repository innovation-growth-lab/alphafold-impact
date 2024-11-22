"""
This is a boilerplate pipeline 'data_collection_pdb'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import fetch_pdb_uniprot_map, fetch_uniprot_data


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    return pipeline(
        [
            node(
                func=fetch_pdb_uniprot_map,
                inputs={
                    "pdb_data": "pdb.entries.raw",
                    },
                outputs="uniprot.pdb.map",
                name="fetch_pdb_uniprot_map",
            ),
            node(
                func=fetch_uniprot_data,
                inputs={
                    "uniprot_ids": "uniprot.pdb.map",
                    },
                outputs="uniprot.entries.raw",
                name="fetch_uniprot_details",
            ),
        ],
        tags=["data_collection_uniprot"],
    )
