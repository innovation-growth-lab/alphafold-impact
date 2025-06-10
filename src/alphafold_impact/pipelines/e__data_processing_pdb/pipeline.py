"""
This is a boilerplate pipeline 'data_processing_pdb'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    collect_pdb_details,
    merge_uniprot_data,
    process_similarity_data,
    aggregate_foldseek_to_pdb_level
)


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
                aggregate_foldseek_to_pdb_level,
                inputs={"similarity_chunks": "foldseek.pdb_similarities.raw"},
                outputs="foldseek.pdb_similarities.intermediate",
                name="aggregate_foldseek_to_pdb_level",
            ),
            node(
                process_similarity_data,
                inputs={
                    "pdb_df": "pdb.entries.intermediate",
                    "foldseek_df": "foldseek.pdb_similarities.intermediate",
                    "rcsb_df": "pdb.structure_matches.raw",
                },
                outputs="pdb.enhanced_entries.intermediate",
                name="process_pdb_similarity_data",
            ),
            node(
                merge_uniprot_data,
                inputs={
                    "pdb_df": "pdb.enhanced_entries.intermediate",
                    "uniprot_df": "uniprot.entries.raw",
                },
                outputs="pdb.entries.primary",
                name="merge_uniprot_data",
            ),
        ],
        tags=["data_processing_pdb"],
    )
