from tests.conftest import config_loader, project_context
import pandas as pd
import pytest
import pdb

from alphafold_impact.pipelines.data_collection_gtr.nodes import (
    fetch_gtr_data,
    GtRDataPreprocessor,
)
from alphafold_impact.settings import DYNAMIC_PIPELINES_MAPPING

# ENDPOINTS = [x[0] for x in DYNAMIC_PIPELINES_MAPPING["gtr"]]


@pytest.fixture
def params(project_context):
    """Get the parameters for the GtR API."""
    return project_context.config_loader["parameters"]["gtr"]


class TestNodes:
    """Test the data collection GTR pipeline nodes."""

    @pytest.mark.parametrize("endpoint, label", DYNAMIC_PIPELINES_MAPPING["gtr"])
    def test_pipeline(self, params, endpoint, label):
        """Test that data is fetched from the GtR API."""

        data = fetch_gtr_data(
            parameters=params[endpoint]["param_requests"],
            endpoint=label,
            test=True,
        )

        self._response_is_list(data)
        self._response_object_is_dict(data)
        self._response_dict_has_keys(data)

        data = pd.DataFrame(data)
        preprocess = GtRDataPreprocessor()
        data = preprocess.methods[label](data)

        self._dataframe_is_returned(data)
        self._dataframe_has_no_links_column(data)

    def _response_is_list(self, data):
        assert isinstance(data, list)

    def _response_object_is_dict(self, data):
        assert isinstance(data[0], dict)

    def _response_dict_has_keys(self, data):
        assert "id" in data[0].keys()
        assert "page_fetched_from" in data[0].keys()

    def _dataframe_is_returned(self, data):
        assert isinstance(data, pd.DataFrame)

    def _dataframe_has_no_links_column(self, data):
        assert "links" not in data.columns
