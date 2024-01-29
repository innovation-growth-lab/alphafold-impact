"""
This is a boilerplate pipeline 'data_collection_lens'
generated using Kedro 0.19.1
"""

from kedro.pipeline import Pipeline, pipeline, node
from alphafold_impact import settings
from .nodes import create_request_form, fetch_lens_data


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=unused-argument
    """Create the data collection pipeline for Lens patents.

    Returns:
        Pipeline: The data collection pipeline for Lens patents.
    """
    template_pipeline = pipeline(
        [
            node(
                func=create_request_form,
                inputs=[
                    "params:api.config",
                    "params:month",
                    "params:year",
                    "params:jurisdiction",
                ],
                outputs=["request_body", "headers"],
                name="create_request_form",
            ),
            node(
                func=fetch_lens_data,
                inputs=["params:api.config", "request_body", "headers"],
                outputs="raw",
                name="fetch_lens_data",
            ),
        ]
    )

    pipelines = []
    for jurisdiction, month, year in settings.DYNAMIC_PIPELINES_MAPPING["lens"]:
        pipelines.append(
            pipeline(
                template_pipeline,
                parameters={
                    "params:api.config": "lens.data_collection.api",
                    "params:month": f"lens.data_collection.months.{month}",
                    "params:year": f"lens.data_collection.years.{year}",
                    "params:jurisdiction": f"lens.data_collection.jurisdictions.{jurisdiction}",
                },
                namespace=f"lens.data_collection.{jurisdiction}.{year}.{month}",
                tags=["lens"],
            )
        )
    return sum(pipelines)
