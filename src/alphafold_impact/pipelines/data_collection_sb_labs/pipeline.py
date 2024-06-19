"""
This is a boilerplate pipeline 'data_collection_sb_labs'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    get_candidate_authors,
    calculate_lab_determinants,
    create_candidates_map,
    assign_lab_label,
    get_institution_info
)
from ..data_collection_labs.nodes import (  # pylint: disable=E0402
    fetch_author_publications,
    get_publications_from_labs,
)
from ..data_processing_labs.nodes import combine_lab_results  # pylint: disable=E0402


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    sb_candidate_authors_pipeline = pipeline(
        [
            node(
                func=get_candidate_authors,
                inputs={
                    "alphafold_data": "lab.data_collection.level.0.raw",
                    "baseline_data": "oa.data_collection.subfield.structural_biology.depth.0.intermediate",
                    "seed_data": "oa.data_processing.subfield.structural_biology.primary",
                    "ct_data": "chains.seed_technologies.intermediate",
                },
                outputs=["authors", "alphafold_authors", "ct_ai_authors", "ct_noai_authors", "ct_authors", "other_authors"],
                tags=["process_sb_candidates", "candidate_authors.get_map"],
            ),
            node(
                func=fetch_author_publications,
                inputs={
                    "author_ids": "authors",
                    "from_publication_date": "params:labs.data_collection.from_date",
                    "api_config": "params:labs.data_collection.api",
                },
                outputs="sb_lab.data_collection.candidates.publications.intermediate",
            ),
            node(
                func=calculate_lab_determinants,
                inputs={
                    "dict_loader": "sb_lab.data_collection.candidates.publications.intermediate",
                    "alphafold_authors": "alphafold_authors",
                    "ct_ai_authors": "ct_ai_authors",
                    "ct_noai_authors": "ct_noai_authors",
                    "other_authors": "other_authors",
                },
                outputs="sb_lab.data_collection.candidates.scores.intermediate",
                tags="process_sb_candidates",
            ),
            node(
                func=combine_lab_results,
                inputs={
                    "dict_loader": "sb_lab.data_collection.candidates.scores.intermediate"
                },
                outputs="sb_lab.data_collection.candidates.scores.primary",
                tags=["process_sb_candidates", "combine_sb_scores"],
            ),
            node(
                func=assign_lab_label,
                inputs={
                    "candidate_data": "sb_lab.data_collection.candidates.scores.primary",
                },
                outputs="sb_lab.data_collection.assignment.primary",
                tags=["process_sb_candidates", "assign_lab_label"],
            ),
            node(
                func=get_publications_from_labs,
                inputs={
                    "data": "sb_lab.data_collection.assignment.primary",
                    "from_publication_date": "params:labs.data_collection.from_date",
                    "api_config": "params:labs.data_collection.api",
                },
                outputs="sb_lab.data_collection.publications.raw",
                tags=["get_publications_from_labs", "get_publications_from_sb_labs"],
            ),
        ],
        tags="candidate_authors.sb_labs",
    )

    create_map_pipeline = pipeline(
        [
            node(
                func=create_candidates_map,
                inputs={
                    "alphafold_authors": "alphafold_authors",
                    "ct_ai_authors": "ct_ai_authors",
                    "ct_noai_authors": "ct_noai_authors",
                    "other_authors": "other_authors",
                },
                outputs="sb_lab.data_collection.candidates.map",
            ),
        ],
        tags="candidate_authors.get_map",
    )

    collect_institution_info = pipeline(
        [
            node(
                func=get_institution_info,
                inputs={
                    "author_ids": "sb_lab.data_collection.assignment.primary",
                },
                outputs="sb_lab.data_collection.institution_info.primary",
            ),
        ],
        tags="get_institution_info",
    )
    return sb_candidate_authors_pipeline + create_map_pipeline + collect_institution_info
