"""
This is a boilerplate pipeline 'data_collection_ecr'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    get_unique_authors,
    fetch_candidate_ecr_status,
    fetch_ecr_outputs,
    merge_ecr_data,
)
from ..f_data_collection_foundational_labs.nodes import (  # pylint: disable=E0402
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
                    "publications_data": "publications.data.outputs",
                },
                outputs="ecr.candidate_authors.raw",
                name="get_unique_authors",
            ),
            node(
                func=fetch_candidate_ecr_status,
                inputs={
                    "authors": "ecr.candidate_authors.raw",
                    "from_publication_date": "params:ecr.data_collection.from_author_date",
                    "to_publication_date": "params:ecr.data_collection.to_author_date",
                    "api_config": "params:ecr.data_collection.api",
                },
                outputs="ecr.authors.raw",
                name="fetch_author_ecr_status",
            ),
            node(
                func=fetch_ecr_outputs,
                inputs={
                    "authors": "ecr.authors.raw",
                    "from_publication_date": "params:ecr.data_collection.from_author_date",
                    "api_config": "params:ecr.data_collection.api",
                },
                tags="debug",
                outputs="ecr.publications.raw",
                name="fetch_ecr_outputs",
            ),
            node(
                func=get_institution_info,
                inputs={
                    "author_ids": "ecr.candidate_authors.raw",
                },
                outputs="ecr.institutions.raw",
                name="get_ecr_institution_info",
            ),
            node(
                func=merge_ecr_data,
                inputs={
                    "data_loaders": "ecr.publications.raw",
                    "candidate_authors": "ecr.candidate_authors.raw",
                    "institutions": "ecr.institutions.raw",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "patents_data": "lens.data_processing.primary",
                    "pdb_submissions": "pdb.entries.intermediate",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="ecr.publications.primary",
                name="merge_ecr_data",
            ),
        ]
    )
