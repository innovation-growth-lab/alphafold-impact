"""This module defines the data collection pipeline for the labs data collection.

The pipeline is used to collect papers from the OpenAIRE API, extract candidate
authors from the papers, and fetch publications for the candidate authors.

To run this pipeline, use the following command:

    $ kedro run --pipeline labs_data_collection
"""

from kedro.pipeline import Pipeline, pipeline, node
from alphafold_impact import settings
from .nodes import (
    get_candidate_authors,
    merge_candidate_authors,
    fetch_author_publications,
    get_pi_id

)


def create_pipeline(  # pylint: disable=unused-argument&missing-function-docstring
    **kwargs,
) -> Pipeline:
    template_candidate_authors_pipeline = pipeline(
        [
            node(
                func=get_candidate_authors,
                inputs={"data": "raw", "baseline": "baseline"},
                outputs="authors",
            )
        ]
    )
    subfield_candidate_authors_pipelines = [
        pipeline(
            template_candidate_authors_pipeline,
            inputs={
                "raw": f"lab.data_collection.subfield.{subfield}.raw",
                "baseline": f"params:labs.data_collection.filter.{subfield}",
            },
            outputs={
                "authors": f"lab.data_collection.candidates.{subfield}.intermediate"
            },
            tags=[
                f"candidate_authors.{subfield}",
                "candidate_authors",
                "candidate_authors.subfields",
            ],
        )
        for subfield in settings.DYNAMIC_PIPELINES_MAPPING["oa"]["subfields"]
    ]
    level_candidate_authors_pipelines = [
        pipeline(
            template_candidate_authors_pipeline,
            inputs={
                "raw": f"lab.data_collection.level.{level}.raw",
                "baseline": f"params:labs.data_collection.filter.{level}",
            },
            outputs={"authors": f"lab.data_collection.candidates.{level}.intermediate"},
            tags=[
                f"candidate_authors.{level}",
                "candidate_authors",
                "candidate_authors.levels",
            ],
        )
        for level in settings.DYNAMIC_PIPELINES_MAPPING["depth_levels"]
    ]

    merge_candidate_authors_pipeline = pipeline(
        [
            node(
                func=merge_candidate_authors,
                inputs={"data": "lab.data_collection.candidates.intermediate"},
                outputs="lab.data_collection.candidates",
            )
        ],
        tags=["merge_candidate_authors", "fetch_author_publications"],
    )

    fetch_author_publications_pipeline = pipeline(
        [
            node(
                func=fetch_author_publications,
                inputs={
                    "author_ids": "lab.data_collection.candidates",
                    "from_publication_date": "params:labs.data_collection.from_date",
                    "api_config": "params:labs.data_collection.api",
                },
                outputs="lab.data_collection.candidates.publications.intermediate",
            )
        ],
        tags=["fetch_author_publications"],
    )

    fetch_ground_truth_author_ids = pipeline(
        [
            node(
                func=get_pi_id,
                inputs={"data": "lab.data_collection.candidates.nih"},
                outputs="lab.data_collection.candidates.nih.intermediate"
            )
        ],
        tags=["fetch_ground_truth_author_ids"],
    )

    return (
        sum(subfield_candidate_authors_pipelines)
        + sum(level_candidate_authors_pipelines)
        + merge_candidate_authors_pipeline
        + fetch_author_publications_pipeline
        + fetch_ground_truth_author_ids
    )
