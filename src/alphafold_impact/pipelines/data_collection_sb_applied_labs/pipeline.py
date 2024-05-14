"""
This is a boilerplate pipeline 'data_collection_other_labs'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (  # pylint: disable=E0402
    get_candidate_authors,
    calculate_lab_determinants,
    create_candidates_map,
    assign_lab_label,
)
from ..data_collection_labs.nodes import (  # pylint: disable=E0402
    fetch_author_publications,
    get_publications_from_labs,
)
from ..data_processing_labs.nodes import combine_lab_results  # pylint: disable=E0402


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    other_candidate_authors_pipeline = pipeline(
        [
            node(
                func=get_candidate_authors,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.primary",
                    "ct_data": "oa.data_processing.structural_biology.depth.ct.intermediate",
                    "other_data": "oa.data_processing.structural_biology.depth.other.intermediate",
                },
                outputs=["applied_authors", "applied_alphafold_authors", "applied_ct_authors", "applied_other_authors"],
                tags=["process_other_candidates", "other_candidate_authors.get_map"],
            ),
            node(
                func=fetch_author_publications,
                inputs={
                    "author_ids": "applied_authors",
                    "from_publication_date": "params:labs.data_collection.from_date",
                    "api_config": "params:labs.data_collection.api",
                },
                outputs="other_lab.data_collection.candidates.publications.intermediate",
            ),
            node(
                func=calculate_lab_determinants,
                inputs={
                    "dict_loader": "other_lab.data_collection.candidates.publications.intermediate",
                    "alphafold_authors": "applied_alphafold_authors",
                    "ct_authors": "applied_ct_authors",
                    "other_authors": "applied_other_authors",
                },
                outputs="other_lab.data_collection.candidates.scores.intermediate",
                tags="process_other_candidates",
            ),
            node(
                func=combine_lab_results,
                inputs={
                    "dict_loader": "other_lab.data_collection.candidates.scores.intermediate"
                },
                outputs="other_lab.data_collection.candidates.scores.primary",
                tags=["process_other_candidates", "combine_other_scores"],
            ),
            node(
                func=assign_lab_label,
                inputs={
                    "candidate_data": "other_lab.data_collection.candidates.scores.primary",
                },
                outputs="other_lab.data_collection.assignment.primary",
                tags=["process_other_candidates", "assign_lab_label"],
            ),
            node(
                func=get_publications_from_labs,
                inputs={
                    "data": "other_lab.data_collection.assignment.primary",
                    "from_publication_date": "params:labs.data_collection.from_date",
                    "api_config": "params:labs.data_collection.api",
                },
                outputs="other_lab.data_collection.publications.raw",
                tags=["get_publications_from_labs"],
            ),
        ],
        tags="candidate_authors.other_labs",
    )

    create_map_pipeline = pipeline(
        [
            node(
                func=create_candidates_map,
                inputs={
                    "alphafold_authors": "applied_alphafold_authors",
                    "ct_authors": "applied_ct_authors",
                    "other_authors": "applied_other_authors",
                },
                outputs="other_lab.data_collection.candidates.map",
            ),
        ],
        tags="other_candidate_authors.get_map",
    )
    return other_candidate_authors_pipeline + create_map_pipeline
