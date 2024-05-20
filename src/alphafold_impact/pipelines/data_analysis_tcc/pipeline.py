"""Pipeline for data anaysis.

This pipeline combines OpenAlex papers with iCite/Pubmed clinical articles.
To run this pipeline, use the following command:

    $ kedro run --pipeline data_analysis_tcc
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import make_tcc_datasets


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=make_tcc_datasets,
                inputs={
                    "icite_data": "nih.data_processing.icite.intermediate",
                    "af_oa_data": "oa.data_processing.depth.intermediate",
                    "biochemistry_oa_data": "oa.data_processing.subfield.biochemistry.primary",
                    "bioinformatics_oa_data": "oa.data_processing.subfield.bioinformatics.primary",
                    "protein_design_oa_data": "oa.data_processing.subfield.protein_design.primary",
                    "structural_biology_oa_data": "oa.data_processing.subfield.structural_biology.primary",
                    "af_label": "params:tcc.af_label",
                    "biochemistry_label": "params:tcc.biochemistry_label",
                    "bioinformatics_label": "params:tcc.bioinformatics_label",
                    "protein_design_label": "params:tcc.protein_design_label",
                    "structural_biology_label": "params:tcc.structural_biology_label",
                },
                outputs="tcc.data_analysis",
            ),
        ],
    )
