"""
This is a boilerplate pipeline 'data_collection_applied_labs'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (  # pylint: disable=E0402
    get_candidate_authors,
    calculate_lab_determinants,
    create_candidates_map,
)
from ..f_data_collection_foundational_labs.nodes import (  # pylint: disable=E0402
    fetch_author_publications,
    get_publications_from_labs,
    get_institution_info,
    combine_lab_results,
    assign_lab_label,
)


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
    **kwargs,
) -> Pipeline:
    other_candidate_authors_pipeline = pipeline(
        [
            node(
                func=get_candidate_authors,
                inputs={
                    "alphafold_data": "oa.data_processing.depth.all.primary",
                    "ct_data": "oa.data_processing.structural_biology.depth.ct.intermediate",
                    "other_data": "oa.data_processing.structural_biology.depth.other.intermediate",
                    "seed_papers": "chains.seed_technologies.intermediate",
                },
                outputs=[
                    "applied_authors",
                    "applied_alphafold_authors",
                    "applied_ct_authors",
                    "applied_ct_ai_authors",
                    "applied_ct_noai_authors",
                    "applied_other_authors",
                ],
                name="get_candidate_applied_authors",
            ),
            node(
                func=fetch_author_publications,
                inputs={
                    "author_ids": "applied_authors",
                    "from_publication_date": "params:labs.data_collection.from_date",
                    "api_config": "params:labs.data_collection.api",
                },
                outputs="applied_lab.data_collection.candidates.publications.intermediate",
                name="fetch_candidate_applied_publications",
            ),
            node(
                func=calculate_lab_determinants,
                inputs={
                    "dict_loader": "applied_lab.data_collection.candidates.publications.intermediate", # pylint: disable=line-too-long
                    "alphafold_authors": "applied_alphafold_authors",
                    "ct_ai_authors": "applied_ct_ai_authors",
                    "ct_noai_authors": "applied_ct_noai_authors",
                    "other_authors": "applied_other_authors",
                },
                outputs="applied_lab.data_collection.candidates.scores.intermediate",
                name="calculate_applied_lab_determinants",
            ),
            node(
                func=combine_lab_results,
                inputs={
                    "dict_loader": "applied_lab.data_collection.candidates.scores.intermediate"
                },
                outputs="applied_lab.data_collection.candidates.scores.primary",
                name="combine_applied_lab_results",
            ),
            node(
                func=assign_lab_label,
                inputs={
                    "candidate_data": "applied_lab.data_collection.candidates.scores.primary",
                    "quantile_val": "params:labs.quantile.applied",
                },
                outputs="applied_lab.data_collection.assignment.primary",
                name="assign_applied_lab_label",
            ),
            node(
                func=get_publications_from_labs,
                inputs={
                    "data": "applied_lab.data_collection.assignment.primary",
                    "from_publication_date": "params:labs.data_collection.from_date",
                    "api_config": "params:labs.data_collection.api",
                },
                outputs="applied_lab.data_collection.publications.raw",
                name="get_publications_from_applied_labs",
            ),
        ],
        tags=["data_collection_applied_labs"],
    )

    create_map_pipeline = pipeline(
        [
            node(
                func=create_candidates_map,
                inputs={
                    "alphafold_authors": "applied_alphafold_authors",
                    "ct_ai_authors": "applied_ct_ai_authors",
                    "ct_noai_authors": "applied_ct_noai_authors",
                    "other_authors": "applied_other_authors",
                },
                outputs="applied_lab.data_collection.candidates.map",
                name="create_applied_candidates_map",
            ),
        ],
        tags=["data_collection_applied_labs"],
    )

    collect_institution_info = pipeline(
        [
            node(
                func=get_institution_info,
                inputs={
                    "author_ids": "applied_lab.data_collection.assignment.primary",
                },
                outputs="applied_lab.data_collection.institution_info.primary",
                name="get_applied_institution_info",
            ),
        ],
        tags=["data_collection_applied_labs"],
    )

    return (
        other_candidate_authors_pipeline
        + create_map_pipeline
        + collect_institution_info
    )
