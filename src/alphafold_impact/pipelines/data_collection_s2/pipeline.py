"""
This is a boilerplate pipeline 'data_collection_s2'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, node, pipeline
from alphafold_impact import settings
from .nodes import (
    load_alphafold_citation_ids,
    get_alphafold_citation_details,
    fetch_citation_strength,
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613

    references_to_alphafold = pipeline(
        [
            node(
                func=load_alphafold_citation_ids,
                inputs={"input_loaders": "cites.input", "work_id": "params:af_paper"},
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
                    "perpage": "params:api.perpage",
                },
                outputs=["level.0", "level.0.log.ids.out"],
            ),
        ],
        tags=["s2.level.0", "strength"],
        namespace="s2.data_collection.strength",
    )

    citation_strength = pipeline(
        [
            node(
                func=fetch_citation_strength,
                inputs={
                    "parent_data": "strength.parent.level",
                    "logs": "log.ids.in",
                    "base_url": "params:api.base_url",
                    "direction": "params:api.cites",
                    "fields": "params:api.fields",
                    "perpage": "params:api.perpage",
                },
                outputs=["raw", "log.ids.out"],
            ),
        ],
    )

    level_pipelines = [
        pipeline(
            citation_strength,
            inputs={
                "strength.parent.level": f"s2.data_collection.strength.level.{level}",
                "log.ids.in": f"s2.data_collection.strength.level.{level}.log.ids.in",
            },
            parameters={
                "api.base_url": "params:s2.data_collection.strength.api.base_url",
                "api.cites": "params:s2.data_collection.strength.api.cites",
                "api.fields": "params:s2.data_collection.strength.api.fields",
                "api.perpage": "params:s2.data_collection.strength.api.perpage",
            },
            namespace=f"s2.data_collection.strength.level.{level+1}",
            tags=["s2.levels", "strength"],
        )
        for level in settings.DYNAMIC_PIPELINES_MAPPING["s2"]
    ]

    return references_to_alphafold + sum(level_pipelines)
