"""
This is a boilerplate pipeline 'data_collection_ecr'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import get_unique_authors, fetch_author_outputs
from ..data_collection_sb_labs.nodes import (  # pylint: disable=E0402
    get_institution_info,
)


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    return pipeline(
        [
            node(
                func=get_unique_authors,
                inputs={
                    "alphafold_data": "lab.data_collection.level.0.raw",
                    "baseline_data": "oa.data_collection.subfield.structural_biology.depth.0.intermediate",
                    "seed_data": "oa.data_processing.subfield.structural_biology.primary",
                    "ct_data": "chains.seed_technologies.intermediate",
                },
                outputs="ecr.authors.raw",
                tags=["process_ecr_ids"],
            ),
            node(
                func=fetch_author_outputs,
                inputs={
                    "author_ids": "ecr.authors.raw",
                    "from_publication_date": "params:labs.data_collection.from_author_date",
                    "api_config": "params:labs.data_collection.api",
                },
                outputs="ecr.authors.publications.raw",
                tags=["collect_ecr_pubs"],
            ),
            node(
                func=get_institution_info,
                inputs={
                    "author_ids": "ecr.authors.raw",
                },
                outputs="ecr.authors.institutions.raw",
                tags="ecr_institutions",
            ),
        ]
    )
