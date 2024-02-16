"""
This module contains functions for fetching NIH funding opportunities data.

Functions:
    collect_nih_funding_opportunities: Downloads and saves NIH
        funding opportunities to S3.
"""

import logging
import time
from typing import Dict, Union
import requests
from requests.adapters import HTTPAdapter, Retry
import pandas as pd

logger = logging.getLogger(__name__)


def _get_params(end_date: str, perpage: int) -> Dict[str, Union[int, str]]:
    """
    Returns the parameters for the NIH grants API requests with
    specified end date and records per page.

    Args:
        end_date (str): The end date for the date range filter in
            the format 'DDMMYYYY'.
        perpage (int): Number of records per page.

    Returns:
        Dict[str, Union[int, str]]: A dictionary of parameters for the API request.
    """
    return {
        "perpage": perpage,
        "sort": "reldate:desc",
        "type": "active,expired,activenosis,expirednosis",
        "parentic": "all",
        "primaryic": "all",
        "activitycodes": "all",
        "doctype": "all",
        "parentfoa": "all",
        "daterange": f"01021991-{end_date}",
        "clinicaltrials": "all",
        "fields": "all",
        "spons": "true",
        "query": "",
    }


def _make_api_request(
    url: str, params: Dict[str, Union[int, str]]
) -> Union[Dict[str, any], None]:
    """
    Creates a session, makes an API request, and returns the
    JSON response.

    Args:
        url (str): The API endpoint URL.
        params (Dict[str, Union[int, str]]) : Parameters for the API request.

    Returns:
        Union[Dict[str, any], None]: The JSON response as a dictionary,
            or None if the request fails.
    """
    session = requests.Session()
    retries = Retry(
        total=5, backoff_factor=0.3, status_forcelist=[429, 500, 502, 503, 504]
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))

    try:
        response = session.get(url, params=params, timeout=10)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        logger.error("Request failed: %s", e)
        return None


def _fetch_total_results(api_url: str, end_date: str) -> int:
    """
    Fetches the total number of NIH grant opportunities up to a
    specified end date.

    Args:
        api_url (str): API URL to collect NIH grant opportunities.
        end_date (str): The end date for the date range filter in the
            format 'DDMMYYYY'.

    Returns:
        int: Total number of results as an integer.
    """
    params = _get_params(end_date, 1)
    params["from"] = 0

    if json_response := _make_api_request(api_url, params):
        return json_response["data"]["hits"]["total"]
    logger.error("Failed to fetch total results count for end date %s", end_date)
    return 0


def _fetch_nih_funding_opportunities(
    api_url: str, end_date: str, total_results: int, perpage: int = 500
) -> pd.DataFrame:
    """
    Fetches all NIH grant opportunities data up to a specified end date.

    Args:
        api_url (str): API URL to collect NIH grant opportunities.
        end_date (str): The end date for the date range filter in the
            format 'DDMMYYYY'.
        total_results (int): Total number of NIH grant opportunities.
        perpage (int): Number of records per page. Defaults to 100.

    Returns:
        pd.DataFrame: A pandas DataFrame containing all NIH grant opportunities.
    """
    params = _get_params(end_date, perpage)
    results = []
    for start in range(0, total_results, perpage):
        params["from"] = start
        if json_response := _make_api_request(api_url, params):
            data = pd.json_normalize(json_response["data"]["hits"]["hits"])
            results.append(data)
            logger.info(
                "Collecting NIH funding opportunities from %d/%d", start, total_results
            )
        else:
            logger.error("Failed to fetch data for start=%d", start)

        time.sleep(1)

    return pd.concat(results, ignore_index=True)


def collect_nih_funding_opportunities(api_url: str, end_date: str) -> pd.DataFrame:
    """
    Fetches and processes NIH grant opportunities data up to today's date.

    Args:
        api_url (str): API URL to collect NIH grant opportunities.
        end_date (str): Latest date to collect funding opportunities from
            in the format ddmmyyyy.

    Returns:
        pd.DataFrame: A pandas DataFrame containing all NIH grant opportunities up
            to today's date.
    """
    try:
        total_results = _fetch_total_results(api_url, end_date)
        return _fetch_nih_funding_opportunities(api_url, end_date, total_results)
    except Exception as e:
        logger.error("An error occurred: %s", e)
        return pd.DataFrame()
