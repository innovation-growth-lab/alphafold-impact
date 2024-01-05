"""
Pipeline for collecting incoming (cites) or outgoing (references) citations for a given work.

The main pipeline is a template pipeline that is used to create two pipelines, using the 
`modular_pipeline` function. The two pipelines are `incoming_pipeline` and
`outgoing_pipeline`, which are used to collect incoming and outgoing citations, respectively.

Note that the output folders are also dynamically named based on the direction of the
citations.

To run this pipeline, use the following command:
    
        $ kedro run --pipeline data_collection_oa

To run this pipeline with a specific work ID, use the following command:

        $ kedro run --pipeline data_collection_oa --params=get.work_id=<work_id> 
        (ie. --params=get.work_id=W2741809807)

To run this pipeline with a specific direction, use the following command:

        $ kedro run --pipeline data_collection_oa --tags=<direction> (ie. --tags=incoming)

The pipeline also includes a `downstream_impact_pipeline` that is used to collect the
incoming citations for the papers that cite the given work ID. This pipeline is used to
calculate the downstream impact of a given work by collecting data for all subsequent
work that cites the work_id.

To run this pipeline, use the following command:
    
            $ kedro run --pipeline data_collection_oa --tags=downstream_impact
"""

from kedro.pipeline import Pipeline, node, pipeline
from kedro.pipeline.modular_pipeline import pipeline as mpl

from .nodes import collect_citation_papers, load_work_ids

template_pipeline = pipeline(
    [
        node(
            func=collect_citation_papers,
            inputs={
                "mailto": "params:api.mailto",
                "perpage": "params:api.perpage",
                "work_ids": "params:get.work_id",
                "direction": "params:direction.template",
            },
            outputs="template",
        )
    ]
)

downstream_impact_pipeline = pipeline(
    [
        node(
            func=load_work_ids,
            inputs={
                "work_id": "params:get.work_id",
                "dataset": "works_incoming_citations",
            },
            outputs="work_ids",
        ),
        node(
            func=collect_citation_papers,
            inputs={
                "mailto": "params:api.mailto",
                "perpage": "params:api.perpage",
                "work_ids": "work_ids",
                "direction": "params:direction.incoming",
            },
            outputs="downstream_incoming_citations",
        ),
    ],
    tags="downstream_impact",
)

incoming_pipeline = mpl(
    pipe=template_pipeline,
    parameters={"params:direction.template": "params:direction.incoming"},
    outputs={"template": "works_outgoing_citations"},
    tags="incoming",
)

outgoing_pipeline = mpl(
    pipe=template_pipeline,
    parameters={"params:direction.template": "params:direction.outgoing"},
    outputs={"template": "works_incoming_citations"},
    tags="outgoing",
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    return incoming_pipeline + outgoing_pipeline + downstream_impact_pipeline
