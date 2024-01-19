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
    load_work_ids,
    preprocess_publication_doi,
    create_list_doi_inputs,
    load_referenced_work_ids,
    retrieve_oa_works_for_concepts_and_years
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    template_pipeline = pipeline(
        [
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:mailto",
                    "perpage": "params:perpage",
                    "work_ids": "params:work_id",
                    "filter_criteria": "params:direction",
                },
                outputs="raw",
            )
        ]
    )

    pipelines = []
    for prefix, direction in settings.DYNAMIC_PIPELINES_MAPPING["oa"]:
        pipelines.append(
            pipeline(
                template_pipeline,
                parameters={
                    "params:direction": f"params:{prefix}.direction.{direction}",
                    "params:mailto": f"params:{prefix}.api.mailto",
                    "params:perpage": f"params:{prefix}.api.perpage",
                    "params:work_id": f"params:{prefix}.get.work_id",
                },
                outputs={"raw": f"{prefix}.{direction}"},
                namespace=f"{prefix}.{direction}",
                tags=[direction, "oa"],
            )
        )

    downstream_impact_pipeline = pipeline(
        [
            node(
                func=load_work_ids,
                inputs={
                    "work_id": "params:get.work_id",
                    "dataset": "cites",
                },
                outputs="work_ids",
            ),
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:api.mailto",
                    "perpage": "params:api.perpage",
                    "work_ids": "work_ids",
                    "filter_criteria": "params:direction.cites",
                },
                outputs="downstream.cites",
            ),
        ],
        namespace="oa.data_collection",
        tags="downstream_impact",
    )


    works_for_concepts_and_years = pipeline(
        [
            node(
                func=retrieve_oa_works_for_concepts_and_years,
                inputs=[
                    "params:test_concept_ids",
                    "params:test_publication_years",
                ],
                outputs="oa_raw_works_for_concepts_and_years",
            ),
        ],
        tags="works_for_concepts_and_years",
    )

    gtr_collection_pipeline = pipeline(
        [
            node(
                func=preprocess_publication_doi,
                inputs="publications",
                outputs="preproc",
                tags=["gtr.citations", "gtr.publications"],
            ),
            node(
                func=create_list_doi_inputs,
                inputs="preproc",
                outputs="doi_list",
                tags=["gtr.citations", "gtr.publications"],
            ),
            node(
                func=collect_papers,
                inputs={
                    "mailto": "params:api.mailto",
                    "perpage": "params:api.perpage",
                    "work_ids": "doi_list",
                    "filter_criteria": "params:filter_doi",
                    "group_work_ids": "params:config.group_work_ids",
                    "slice_keys": "params:config.slice_keys",
                    "parallelise": "params:config.parallelise",
                },
                outputs="works",
                tags=["gtr.citations", "gtr.publications"],
            ),
            node(
                func=load_referenced_work_ids,
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
                    "work_ids": "oa_list",
                    "filter_criteria": "params:filter_oa",
                    "group_work_ids": "params:config.group_work_ids",
                    "slice_keys": "params:config.slice_keys",
                    "parallelise": "params:config.parallelise",
                },
                outputs="citations",
                tags=["gtr.citations", "gtr.references"],
            ),
        ],
        namespace="oa.data_collection.gtr",
    )

    return (
        sum(pipelines)  # Base pipelines
        + downstream_impact_pipeline  # Base pipelines
        + gtr_collection_pipeline  # GtR pipelines
        + works_for_concepts_and_years
    )
