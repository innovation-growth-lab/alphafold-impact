"""
This is a boilerplate pipeline 'data_collection_sb_labs'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import get_candidate_authors, calculate_lab_determinants
from ..data_collection_labs.nodes import fetch_author_publications  # pylint: disable=E0402


def create_pipeline(  # pylint: disable=unused-argument,missing-function-docstring
        **kwargs) -> Pipeline:
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
                outputs=["authors", "alphafold_authors", "ct_authors", "other_authors"],
                tags="process_sb_candidates"
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
                    "ct_authors": "ct_authors",
                    "other_authors": "other_authors",
                },
                outputs="sb_lab.data_collection.candidates.scores.intermediate",
                tags="process_sb_candidates"
            ),

        ],
        tags="candidate_authors.sb_labs",
    )
    return sb_candidate_authors_pipeline
