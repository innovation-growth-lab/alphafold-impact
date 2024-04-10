"""Pipeline for reporting.

This pipeline produces TCC charts and tables.
To run this pipeline, use the following command:

    $ kedro run --pipeline reporting_tcc
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import save_tcc_summary_subfields


def aggregate_datasets(*datasets):
    """Aggregates multiple datasets into a list."""
    return list(datasets)


def create_pipeline(**kwargs) -> Pipeline:
    return Pipeline(
        [
            node(
                func=aggregate_datasets,
                inputs=[
                    "oa.data_processing.depth.intermediate",
                    "oa.data_processing.subfield.biochemistry.primary",
                    "oa.data_processing.subfield.bioinformatics.primary",
                    "oa.data_processing.subfield.protein_design.primary",
                    "oa.data_processing.subfield.structural_biology.primary",
                ],
                outputs="aggregated_oa_datasets",
            ),
            node(
                func=save_tcc_summary_subfields,
                inputs={
                    "tcc_subfield_partitions": "tcc.data_analysis",
                    "oa_datasets": "aggregated_oa_datasets",
                    "partition_order": "params:tcc.partition_order",
                },
                outputs="tcc.reporting.summary_subfields",
            ),
        ],
    )
