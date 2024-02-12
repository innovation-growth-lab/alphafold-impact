"""
This is a boilerplate pipeline 'data_collection_s2'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import (
    load_alphafold_citation_ids,
    get_alphafold_citation_details
)

def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    return pipeline(
        [
            node(
                func=load_alphafold_citation_ids,
                inputs={
                    "input_loaders": "cites.input",
                    "work_id": "params:af_paper"
                },
                outputs="af_citations",
            ),
            node(
                func=get_alphafold_citation_details,
                inputs={
                    "work_ids": "af_citations",
                    "af_doi": "params:af_doi",
                    "base_url": "params:api.base_url",
                    "direction": "params:api.cited_by",
                    "fields": "params:api.fields",
                    "perpage": "params:api.perpage"
                },
                outputs="intent.alphafold.raw",
            )
        ],
        namespace="s2.data_collection"
    )
