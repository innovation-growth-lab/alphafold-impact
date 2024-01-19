# # pylint: skip-file
# """
# This module contains unit tests for the data collection
#     OA pipeline nodes.

# The TestNodes class defines test cases for the pipeline nodes
#     used in the OA data collection process.
# """

import pandas as pd
import pytest

from alphafold_impact.pipelines.data_collection_oa.nodes import (
    collect_papers,
)
from alphafold_impact.settings import DYNAMIC_PIPELINES_MAPPING

expected_keys = [
    "id",
    "doi",
    "display_name",
    "title",
    "publication_date",
    "abstract",
    "authorships",
    "cited_by_count",
    "concepts",
    "keywords",
    "grants",
    "referenced_works",
]


def test_data_is_list(data):
    assert isinstance(data, list), "Data should be a list"


def test_items_are_dicts(data):
    for item in data:
        assert isinstance(item, dict), "Each item in data should be a dictionary"


def test_expected_keys_in_dicts(data):
    for item in data:
        assert set(item.keys()) == set(
            expected_keys
        ), "Each item should have the expected keys"


def test_id_is_not_empty(data):
    for item in data:
        assert item["id"], "id should not be empty"


@pytest.fixture
def params(project_context):
    """Get the parameters for the GtR API."""
    return project_context.config_loader["parameters"]["oa"]["data_collection"]


@pytest.mark.oa
@pytest.mark.parametrize("filter_", DYNAMIC_PIPELINES_MAPPING["oa"])
def test_collect_papers(params, filter_):
    """Test that data is fetched from the OA API and processed correctly.

    Args:
        params (dict): The parameters for the GtR API.
        endpoint (str): The API endpoint to fetch data from.
        label (str): The label specifying the data preprocessing method.

    Raises:
        AssertionError: If any of the assertions fail.
    """

    # Define the data collection pipeline
    loader_dict = collect_papers(
        mailto=params["direction"][filter_]["api"]["mailto"],
        perpage=params["direction"][filter_]["api"]["perpage"],
        work_ids=params["direction"][filter_]["get"]["work_id"],
        filter_criteria=params["direction"][filter_]["filter"],
    )

    # Load the data
    data = list(loader_dict.values())[0]()

    # Assert that the response is a list
    test_data_is_list(data)

    # Assert that each item in the response list is a dictionary
    test_items_are_dicts(data)

    # Assert that each dictionary in the response has the required keys
    test_expected_keys_in_dicts(data)

    # Assert that the id is not empty
    test_id_is_not_empty(data)
