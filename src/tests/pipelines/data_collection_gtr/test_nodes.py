# pylint: skip-file
"""
This module contains unit tests for the data collection 
    GTR pipeline nodes.

The TestNodes class defines test cases for the pipeline nodes 
    used in the GTR data collection process.
"""

import pandas as pd
import pytest

from alphafold_impact.pipelines.data_collection_gtr.nodes import (
    fetch_gtr_data,
    GtRDataPreprocessor,
)
from alphafold_impact.settings import DYNAMIC_PIPELINES_MAPPING


@pytest.fixture
def params(project_context):
    """Get the parameters for the GtR API."""
    return project_context.config_loader["parameters"]["gtr"]


class TestNodes:
    """Test the data collection GTR pipeline nodes."""

    @pytest.mark.parametrize("endpoint, label", DYNAMIC_PIPELINES_MAPPING["gtr"])
    def test_pipeline(self, params, endpoint, label):
        """
        Test that data is fetched from the GtR API and processed correctly.

        Args:
            params (dict): The parameters for the GtR API.
            endpoint (str): The API endpoint to fetch data from.
            label (str): The label specifying the data preprocessing method.

        Raises:
            AssertionError: If any of the assertions fail.
        """
        # Fetch data from the GtR API
        data = fetch_gtr_data(
            parameters=params[endpoint]["param_requests"],
            endpoint=label,
            test=True,
        )

        # Assert that the response is a list
        self._response_is_list(data)
        # Assert that each item in the response list is a dictionary
        self._response_object_is_dict(data)
        # Assert that each dictionary in the response has the required keys
        self._response_dict_has_keys(data)

        # Convert the response data into a pandas DataFrame
        data = pd.DataFrame(data)
        # Preprocess the data based on the specified label
        preprocess = GtRDataPreprocessor()
        data = preprocess.methods[label](data)

        # Assert that the processed data is a pandas DataFrame
        self._dataframe_is_returned(data)
        # Assert that the processed data does not have a "links" column
        self._dataframe_has_no_links_column(data)

    def _response_is_list(self, data):
        """
        Assert that the response data is a list.

        Args:
            data (list): The response data to be checked.

        Raises:
            AssertionError: If the data is not a list.
        """
        assert isinstance(data, list)

    def _response_object_is_dict(self, data):
        """
        Assert that each item in the response list is a dictionary.

        Args:
            data (list): The response data to be checked.

        Raises:
            AssertionError: If any item in the data is not a dictionary.
        """
        assert isinstance(data[0], dict)

    def _response_dict_has_keys(self, data):
        """
        Assert that each dictionary in the response has the required keys.

        Args:
            data (list): The response data to be checked.

        Raises:
            AssertionError: If any dictionary in the data is missing the required keys.
        """
        assert "id" in data[0].keys()
        assert "page_fetched_from" in data[0].keys()

    def _dataframe_is_returned(self, data):
        """
        Assert that the processed data is a pandas DataFrame.

        Args:
            data (pd.DataFrame): The processed data to be checked.

        Raises:
            AssertionError: If the data is not a pandas DataFrame.
        """
        assert isinstance(data, pd.DataFrame)

    def _dataframe_has_no_links_column(self, data):
        """
        Assert that the processed data does not have a "links" column.

        Args:
            data (pd.DataFrame): The processed data to be checked.

        Raises:
            AssertionError: If the data has a "links" column.
        """
        assert "links" not in data.columns
