"""
This is a boilerplate pipeline 'data_collection_ecr'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    get_unique_authors,
    fetch_candidate_ecr_status,
    merge_ecr_authors_data,
    fetch_author_outputs,
    merge_author_data,
    aggregate_to_quarterly,
)
from ..f_data_collection_foundational_labs.nodes import (  # pylint: disable=E0402
    get_institution_info,
)


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    basic_pipeline = pipeline(
        [
            # node(
            #     func=get_unique_authors,
            #     inputs={
            #         "publications_data": "publications.data.outputs",
            #     },
            #     outputs="ecr.candidate_authors.raw",
            #     name="get_unique_authors",
            # ),
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
                func=get_institution_info,
                inputs={
                    "author_ids": "ecr.candidate_authors.raw",
                },
                outputs="ecr.institutions.raw",
                name="get_ecr_institution_info",
            ),
            node(
                func=merge_ecr_authors_data,
                inputs={
                    "data_loaders": "ecr.authors.raw",
                },
                outputs="ecr.authors.merged",
                name="merge_ecr_authors_data",
            ),
        ],
        tags=["basic_ecr_collection"],

    )

    ecr_pipeline = pipeline(
        [
            node(
                func=fetch_author_outputs,
                inputs={
                    "authors": "ecr.authors.merged",
                    "ecr": "params:ecr.data_collection.ecr_bool.ecr",
                    "from_publication_date": "params:ecr.data_collection.from_ecr_date",
                    "api_config": "params:ecr.data_collection.api",
                },
                outputs="ecr.publications.raw",
                name="fetch_ecr_outputs",
            ),
            node(
                func=merge_author_data,
                inputs={
                    "data_loaders": "ecr.publications.raw",
                    "candidate_authors": "ecr.candidate_authors.raw",
                    "authors": "ecr.authors.merged",
                    "institutions": "ecr.institutions.raw",
                    "ecr": "params:ecr.data_collection.ecr_bool.ecr",
                    "patents_data": "lens.data_processing.primary",
                    "pdb_submissions": "pdb.entries.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                },
                outputs="ecr.publications.primary",
                name="merge_ecr_data",
            ),
            node(
                func=aggregate_to_quarterly,
                inputs={
                    "data": "ecr.publications.primary",
                },
                outputs="ecr.publications.quarterly",
                name="aggregate_ecr_to_quarterly",
            ),
        ],
        tags=["ecr_pipeline"],
    )

    nonecr_pipeline = pipeline(
        [
            node(
                func=fetch_author_outputs,
                inputs={
                    "authors": "ecr.authors.merged",
                    "ecr": "params:ecr.data_collection.ecr_bool.nonecr",
                    "from_publication_date": "params:ecr.data_collection.from_nonecr_date",
                    "api_config": "params:ecr.data_collection.api",
                },
                outputs="nonecr.publications.raw",
                name="fetch_nonecr_outputs",
            ),
            node(
                func=merge_author_data,
                inputs={
                    "data_loaders": "nonecr.publications.raw",
                    "candidate_authors": "ecr.candidate_authors.raw",
                    "authors": "ecr.authors.merged",
                    "institutions": "ecr.institutions.raw",
                    "ecr": "params:ecr.data_collection.ecr_bool.nonecr",
                    "patents_data": "lens.data_processing.primary",
                    "pdb_submissions": "pdb.entries.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                },
                outputs="nonecr.publications.primary",
                name="merge_nonecr_data",
                # tags="debug",
            ),
            node(
                func=aggregate_to_quarterly,
                inputs={
                    "data": "nonecr.publications.primary",
                },
                outputs="nonecr.publications.quarterly",
                name="aggregate_nonecr_to_quarterly",
                # tags="debug",
            ),
        ],
        tags=["nonecr_pipeline"],
    )

    return (
        basic_pipeline + 
        ecr_pipeline + 
        nonecr_pipeline
    )