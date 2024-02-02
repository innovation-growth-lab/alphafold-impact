"""
Pipeline for collecting Lens patent data. It includes:

    Base:
    - Getting Lens API credentials.
    - Creating a request form for Lens data collection.
    - Fetching Lens data.

To run this pipeline, use the following command:

    $ kedro run --pipeline data_collection_lens

To run a specific pipeline, use a combination of the --namespace and --tags flags:

    $ kedro run --pipeline data_collection_lens --namespace lens.data_collection.<jurisdiction>.<year>.<month>
    $ kedro run --pipeline data_collection_lens --tags lens

"""

from kedro.pipeline import Pipeline, pipeline, node
from alphafold_impact import settings
from .nodes import create_request_form, fetch_lens_data, get_app_credentials


def create_pipeline(**kwargs) -> Pipeline:  # pylint: disable=unused-argument
    """Create the data collection pipeline for Lens patents.

    Returns:
        Pipeline: The data collection pipeline for Lens patents.
    """
    template_pipeline = pipeline(
        [
            node(
                func=get_app_credentials,
                inputs=None,
                outputs="lens_token",
                name="get_app_credentials",
            ),
            node(
                func=create_request_form,
                inputs=[
                    "params:api.config",
                    "params:month",
                    "params:year",
                    "params:jurisdiction",
                    "lens_token"
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
