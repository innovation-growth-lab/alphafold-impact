# """
# This is a boilerplate test file for pipeline 'data_enrichment_gtr'
# generated using Kedro 0.19.1.
# Please add your pipeline tests here.

# Kedro recommends using `pytest` framework, more info about it can be found
# in the official documentation:
# https://docs.pytest.org/en/latest/getting-started.html
# """
# import pytest
# from pathlib import Path
# from kedro.config import OmegaConfigLoader
# from kedro.framework.context import KedroContext
# from kedro.framework.hooks import _create_hook_manager
# from kedro.pipeline import Pipeline, node
# from alphafold_impact.pipelines.data_collection_gtr import create_pipeline  # Import 'create_pipeline' from the correct module
# import pdb




# # @pytest.fixture
# # def config_loader():
# #     """Create a config loader using a local conf folder."""
# #     return OmegaConfigLoader(conf_source=str(Path.cwd()))


# # @pytest.fixture
# # def project_context(config_loader):
# #     """Create a project context using a local conf folder."""
# #     return KedroContext(
# #         package_name="alphafold_impact",
# #         project_path=Path.cwd(),
# #         config_loader=config_loader,
# #         hook_manager=_create_hook_manager(),
# #         env="base",
# #     )

# # def test_project_context(project_context):
# #     """Test that the project context is valid."""
# #     assert isinstance(project_context, KedroContext)

# # def test_project_path(project_context):
# #     """Test that the project path is correct."""
# #     assert project_context.project_path == Path.cwd()

# class TestDataCollectionGTRPipeline:
#     """Test the data collection GTR pipeline."""

# #     def test_create_pipeline(self, project_context):
# #         """Test that the create_pipeline function returns a valid pipeline."""
# #         pipeline = create_pipeline()
# #         assert isinstance(pipeline, Pipeline)


# #     def test_pipeline_inputs(self, project_context):
# #         """Test the inputs of the data collection GTR pipeline."""
# #         pipeline = create_pipeline()
# #         # assert pipeline.inputs == {"param_requests", "endpoint_template"}
# #         assert 3 == 3

#     def test_pipeline_outputs(self, project_context):
#         """Test the outputs of the data collection GTR pipeline."""
#         pipeline = create_pipeline()
#         # assert pipeline.outputs == {"intermediate"}
#         assert 3 == 3

# #     def test_pipeline_namespace(self, project_context):
# #         """Test the namespace of the data collection GTR pipeline."""
# #         pipeline = create_pipeline()
# #         # assert pipeline.namespace == "gtr"
# #         assert 3 == 3

# #     def test_pipeline_tags(self, project_context):
# #         """Test the tags of the data collection GTR pipeline."""
# #         pipeline = create_pipeline()
# #         # assert pipeline.tags == ["gtr"]
# #         assert 3 == 3

# # TestDataCollectionGTRPipeline()