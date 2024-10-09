"""
Pipeline for collecting OpenAlex data. It includes:
    
    Base:
    - Collecting OpenAlex publications given a work_id.
    - Collect references and citations for a work_id.
    - Collect citations for all citations of a work_id.

    GtR:
        - Collecting OpenAlex publications that match the doi values for the Gateway to Research
            publications.
        - Collecting all papers that these publications cite.

To run this pipeline, use the following command:

    $ kedro run --pipeline data_collection_oa

To run a specific pipeline, use a combination of the --namespace and --tags flags:

    $ kedro run --pipeline data_collection_oa --namespace oa.data_collection.gtr 
    $ kedro run --pipeline data_collection_oa --tags downstream_impact



"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import (
    fetch_subfield_baseline,
    fetch_citation_to_specific_depth,
    preprocess_baseline_data,
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613

    fixed_depth_alphafold_pipeline = pipeline(
        [
            node(
                func=fetch_citation_to_specific_depth,
                inputs={
                    "seed_paper": "params:oa.data_collection.depth.get.work_ids",
                    "api_config": "params:oa.data_collection.depth.api",
                    "filter_config": "params:oa.data_collection.depth.filter",
                    "max_depth": "params:oa.data_collection.depth.max_depth",
                },
                outputs="oa.data_collection.triad.depth.level.raw",
                name="fetch_af_citation_to_specific_depth",
            )
        ],
        tags="data_collection_oa",
    )

    structural_biology_pipeline = pipeline(
        [
            node(
                func=fetch_subfield_baseline,
                inputs={
                    "oa_concept_ids": "params:oa.data_collection.subfield.structural_biology.concept_ids",  # pylint: disable=line-too-long
                    "from_publication_date": "params:oa.data_collection.subfield.structural_biology.from_date",  # pylint: disable=line-too-long
                    "api_config": "params:oa.data_collection.subfield.structural_biology.api",
                },
                outputs="oa.data_collection.subfield.structural_biology.raw",
                name="fetch_baseline_sb",
            )
        ],
        tags="data_collection_oa",
    )

    fixed_depth_structural_biology_pipeline = pipeline(
        [
            node(
                func=preprocess_baseline_data,
                inputs={
                    "data": "oa.data_collection.subfield.structural_biology.raw",
                    "alphafold_papers": "params:oa.data_collection.depth.get.work_ids",
                },
                outputs="oa.structural_biology.seed_papers",
                name="preprocess_sb",
            ),
            node(
                func=fetch_citation_to_specific_depth,
                inputs={
                    "seed_paper": "oa.structural_biology.seed_papers",
                    "api_config": "params:oa.data_collection.depth.api",
                    "filter_config": "params:oa.data_collection.depth.filter",
                    "max_depth": "params:sb_labs.data_collection.depth.max_depth",
                },
                outputs="oa.data_collection.subfield.structural_biology.depth.raw",
                name="fetch_sb_citation_to_specific_depth",
            ),
        ],
        tags=["data_collection_oa"],
    )


    return (
        fixed_depth_alphafold_pipeline
        + structural_biology_pipeline  # SB pipeline
        + fixed_depth_structural_biology_pipeline
    )
