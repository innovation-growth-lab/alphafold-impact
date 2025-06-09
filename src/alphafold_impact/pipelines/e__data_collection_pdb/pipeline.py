"""
This is a boilerplate pipeline 'data_collection_pdb'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import fetch_pbd_details, collect_structure_match_scores


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    return pipeline(
        [
            node(
                func=fetch_pbd_details,
                inputs={"config": "params:pdb.config"},
                outputs="pdb.entries.raw",
                name="fetch_pdb_details",
            ),
            node(
                func=collect_structure_match_scores,
                inputs={
                    "pdb_entries": "pdb.entries.raw",
                    "config": "params:pdb.config",
                },
                outputs="pdb.structure_matches.raw",
                name="collect_structure_match_scores",
            ),
        ],
        tags=["data_collection_pdb"],
    )
