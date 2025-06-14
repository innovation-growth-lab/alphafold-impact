"""
This is a boilerplate pipeline 'data_collection_foundational_labs'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    assign_lab_label,
    calculate_lab_determinants,
    combine_lab_results,
    create_candidates_map,
    fetch_author_publications,
    get_candidate_authors,
    get_institution_info,
    get_publications_from_labs,
)


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    foundational_candidate_authors_pipeline = pipeline(
        [
            node(
                func=get_candidate_authors,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.0.intermediate",
                    "baseline_data": "oa.data_collection.subfield.structural_biology.depth.0.intermediate", # pylint: disable=line-too-long
                    "seed_data": "oa.data_processing.subfield.structural_biology.primary",
                    "ct_data": "chains.seed_technologies.intermediate",
                },
                outputs=[
                    "authors",
                    "alphafold_authors",
                    "ct_ai_authors",
                    "ct_pp_authors",
                    "ct_sb_authors",
                    "ct_authors",
                    "other_authors",
                ],
                name="get_candidate_foundational_authors",
            ),
            node(
                func=fetch_author_publications,
                inputs={
                    "author_ids": "authors",
                    "from_publication_date": "params:labs.data_collection.from_date",
                    "api_config": "params:labs.data_collection.api",
                },
                outputs="foundational_lab.data_collection.candidates.publications.intermediate",
                name="fetch_candidate_foundational_publications",
            ),
            node(
                func=calculate_lab_determinants,
                inputs={
                    "dict_loader": "foundational_lab.data_collection.candidates.publications.intermediate",
                    "alphafold_authors": "alphafold_authors",
                    "ct_ai_authors": "ct_ai_authors",
                    "ct_pp_authors": "ct_pp_authors",
                    "ct_sb_authors": "ct_sb_authors",
                    "other_authors": "other_authors",
                },
                outputs="foundational_lab.data_collection.candidates.scores.intermediate",
                name="calculate_foundational_lab_determinants",
            ),
            node(
                func=combine_lab_results,
                inputs={
                    "dict_loader": "foundational_lab.data_collection.candidates.scores.intermediate"
                },
                outputs="foundational_lab.data_collection.candidates.scores.primary",
                name="combine_foundational_lab_results",
            ),
            node(
                func=assign_lab_label,
                inputs={
                    "candidate_data": "foundational_lab.data_collection.candidates.scores.primary",
                    "quantile_val": "params:labs.quantile.foundational",
                },
                outputs="foundational_lab.data_collection.assignment.primary",
                name="assign_foundational_lab_label",
            ),
            node(
                func=get_publications_from_labs,
                inputs={
                    "data": "foundational_lab.data_collection.assignment.primary",
                    "from_publication_date": "params:labs.data_collection.from_date",
                    "api_config": "params:labs.data_collection.api",
                },
                outputs="foundational_lab.data_collection.publications.raw",
                name="get_publications_from_foundational_labs",
            ),
        ],
        tags=["data_collection_foundational_labs"],
    )

    create_map_pipeline = pipeline(
        [
            node(
                func=create_candidates_map,
                inputs={
                    "alphafold_authors": "alphafold_authors",
                    "ct_ai_authors": "ct_ai_authors",
                    "ct_pp_authors": "ct_pp_authors",
                    "ct_sb_authors": "ct_sb_authors",
                    "other_authors": "other_authors",
                },
                outputs="foundational_lab.data_collection.candidates.map",
                name="create_foundational_candidates_map",
            ),
        ],
        tags=["data_collection_foundational_labs"],
    )

    collect_institution_info = pipeline(
        [
            node(
                func=get_institution_info,
                inputs={
                    "author_ids": "foundational_lab.data_collection.assignment.primary",
                },
                outputs="foundational_lab.data_collection.institution_info.primary",
                name="get_foundational_institution_info",
            ),
        ],
        tags=["data_collection_foundational_labs"],
    )
    return (
        foundational_candidate_authors_pipeline
        + create_map_pipeline
        + collect_institution_info
    )
