"""
This is a boilerplate pipeline 'publications_descriptive_translational'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import (
    create_publications_data,
    merge_individual_data,
    update_alphafold_triad,
    define_high_pdb_authors,
    select_regression_columns,
    get_institution_info,
)


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=C0116,W0613
    create_publications_data_pipeline = pipeline(
        [
            node(
                update_alphafold_triad,
                inputs={
                    "data": "oa.chain_labels.id.primary",
                },
                outputs="oa.chain_labels.id.updated",
                name="update_alphafold_triad",
            ),
            node(
                create_publications_data,
                inputs={
                    "data": "oa.chain_labels.id.updated",
                    "source": "params:publications.source.af",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "patents_data": "lens.data_processing.primary",
                    "pdb_submissions": "pdb.entries.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="publications.data.af",
                name="create_publications_data_af",
            ),
            node(
                create_publications_data,
                inputs={
                    "data": "oa.chain_labels.id.ct.primary",
                    "seed_papers": "chains.seed_technologies.intermediate",
                    "source": "params:publications.source.ct",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "patents_data": "lens.data_processing.primary",
                    "pdb_submissions": "pdb.entries.primary",
                    "icite_data": "pubmed.data_processing.icite.intermediate",
                },
                outputs="publications.data.ct",
                name="create_publications_data_ct_sb",
            ),
            node(
                create_publications_data,
                inputs={
                    "data": "oa.chain_labels.id.other.primary",
                    "source": "params:publications.source.other",
                    "mesh_terms": "nih.data_collection.mesh_terms",
                    "patents_data": "lens.data_processing.primary",
                    "pdb_submissions": "pdb.entries.primary",
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
                outputs="publications.data.merged",
                name="merge_publications_data",
                tags="complete",
            ),
            node(
                get_institution_info,
                inputs={"publications": "publications.data.merged"},
                outputs="publications.data.intermediate",
                name="get_last_author_institution_info",
                tags="complete",
            ),
            node(
                define_high_pdb_authors,
                inputs={
                    "data": "publications.data.intermediate",
                    "pdb_submissions": "pdb.entries.primary",
                },
                outputs="publications.data.outputs",
                name="define_high_pdb_authors",
                tags="complete",
            ),
            node(
                select_regression_columns,
                inputs={
                    "data": "publications.data.outputs",
                    "columns": "params:publications.columns_to_drop",
                },
                outputs="publications.regression.inputs",
                name="subset_columns_for_regression",
                tags="complete",
            ),
        ],
        tags=["data_output_publications"],
    )

    return create_publications_data_pipeline
