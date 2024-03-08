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
from alphafold_impact import settings

from .nodes import (
    collect_papers,
    load_oa_ids,
    preprocess_publication_doi,
    create_list_doi_inputs,
    load_referenced_oa_ids,
    retrieve_oa_works_for_concepts_and_years,
    fetch_subfield_baseline,
    fetch_subfield_and_logic,
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    template_citation_pipeline = pipeline(
        [
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:api.mailto",
                    "perpage": "params:api.perpage",
                    "oa_ids": "params:get.work_id",
                    "filter_criteria": "params:filter",
                },
                outputs="raw",
            )
        ]
    )

    baseline_pipelines = [
        pipeline(
            template_citation_pipeline,
            namespace=f"oa.data_collection.direction.{filter_}",
            tags=[filter_, "oa"],
        )
        for filter_ in settings.DYNAMIC_PIPELINES_MAPPING["oa"]["directions"]
    ]

    downstream_impact_pipeline = pipeline(
        [
            node(
                func=load_oa_ids,
                inputs={
                    "oa_id": "params:get.work_id",
                    "dataset": "cites.input",
                },
                outputs="work_ids",
            ),
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:api.mailto",
                    "perpage": "params:api.perpage",
                    "oa_ids": "work_ids",
                    "filter_criteria": "params:filter",
                },
                outputs="cites.intermediate",
            ),
        ],
        namespace="oa.data_collection.downstream",
        tags="downstream",
    )

    works_for_concepts_and_years = pipeline(
        [
            node(
                func=retrieve_oa_works_for_concepts_and_years,
                inputs=[
                    "params:test_concept_ids",
                    "params:test_publication_years",
                ],
                outputs="raw",
            ),
        ],
        namespace="oa.data_collection.works_for_concepts",
        tags="works_for_concepts",
    )

    gtr_collection_pipeline = pipeline(
        [
            node(
                func=preprocess_publication_doi,
                inputs="publications",
                outputs="preproc",
                tags=["gtr.list", "gtr.citations", "gtr.publications"],
            ),
            node(
                func=create_list_doi_inputs,
                inputs="preproc",
                outputs="doi_list",
                tags=["gtr.list", "gtr.citations", "gtr.publications"],
            ),
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:api.mailto",
                    "perpage": "params:api.perpage",
                    "oa_ids": "doi_list",
                    "filter_criteria": "params:filter_doi",
                    "group_oa_ids": "params:config.group_work_ids",
                    "slice_keys": "params:config.slice_keys",
                    "parallelise": "params:config.parallelise",
                },
                outputs="works",
                tags=["gtr.citations", "gtr.publications"],
            ),
            node(
                func=load_referenced_oa_ids,
                inputs="works",
                outputs=[
                    "oa_list",
                    "dict",
                ],
                tags=["gtr.citations", "gtr.references"],
            ),
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:api.mailto",
                    "perpage": "params:api.perpage",
                    "oa_ids": "oa_list",
                    "filter_criteria": "params:filter_oa",
                    "group_oa_ids": "params:config.group_work_ids",
                    "slice_keys": "params:config.slice_keys",
                    "parallelise": "params:config.parallelise",
                },
                outputs="citations",
                tags=["gtr.citations", "gtr.references"],
            ),
        ],
        namespace="oa.data_collection.gtr",
    )

    subfield_baseline_pipeline = pipeline(
        [
            node(
                func=fetch_subfield_baseline,
                inputs={
                    "oa_concept_ids": "params:concept_ids",
                    "from_publication_date": "params:from_date",
                    "api_config": "params:api",
                },
                outputs="raw",
                tags="subfield",
            )
        ]
    )

    subfield_baselines = [
        pipeline(
            subfield_baseline_pipeline,
            namespace=f"oa.data_collection.subfield.{concept}",
            tags=[concept, "subfield"],
        )
        for concept in settings.DYNAMIC_PIPELINES_MAPPING["oa"]["subfields"]
    ]

    subfield_artificial_intelligence_pipeline = pipeline(
        [
            node(
                func=fetch_subfield_and_logic,
                inputs={
                    "oa_main_concept_ids": "params:concept_ids",
                    "oa_and_concept_ids": "params:and_ids",
                    "from_publication_date": "params:from_date",
                    "api_config": "params:api",
                },
                outputs="raw",
            ),
        ],
        tags=["subfield_artificial_intelligence", "subfield"],
        namespace="oa.data_collection.subfield.artificial_intelligence",
    )

    return (
        sum(baseline_pipelines)  # Base pipelines
        + downstream_impact_pipeline  # Papers that cite papers that cite AF
        + works_for_concepts_and_years  # Concept and year pipelines
        + gtr_collection_pipeline  # GtR pipelines
        + sum(subfield_baselines)  # Subfield pipelines
        + subfield_artificial_intelligence_pipeline  # AI pipeline
    )
