"""
This is a boilerplate pipeline 'publications_descriptive_translational'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    create_publications_data,
    merge_individual_data,
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    create_publications_data_pipeline = pipeline(
        [
            node(
                create_publications_data,
                inputs={
                    "data": "oa.chain_labels.id.primary",
                    "source": "params:publications.source.af",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "patents_data": "lens.data_processing.primary",
                    "pdb_submissions": "pdb.entries.intermediate",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="publications.data.af",
                name="create_publications_data_af",
            ),
            node(
                create_publications_data,
                inputs={
                    "data": "oa.chain_labels.id.ct.primary",
                    "source": "params:publications.source.ct",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "patents_data": "lens.data_processing.primary",
                    "pdb_submissions": "pdb.entries.intermediate",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="publications.data.ct",
                name="create_publications_data_ct_sb",
            ),
            node(
                create_publications_data,
                inputs={
                    "data": "oa.chain_labels.id.other.primary",  # "oa.data_processing.structural_biology.depth.other.intermediate",
                    "source": "params:publications.source.other",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "patents_data": "lens.data_processing.primary",
                    "pdb_submissions": "pdb.entries.intermediate",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="publications.data.other",
                name="create_publications_data_other_sb",
            ),
            node(
                merge_individual_data,
                inputs={
                    "data_af": "publications.data.af",
                    "data_ct": "publications.data.ct",
                    "data_other": "publications.data.other",
                },
                outputs="publications.data.outputs",
                name="merge_publications_data",
            ),
        ],
        tags=["rerun", "data_output_publications"],
    )

    return create_publications_data_pipeline