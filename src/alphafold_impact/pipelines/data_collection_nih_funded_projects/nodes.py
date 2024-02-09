"""
This module contains functions for fetching NIH funded projects data.

Functions:
    collect_nih_funded_projects: Downloads and saves NIH funded projects
        to S3 as a partitioned dataset.
"""

import logging
import time
from typing import List, Dict, Tuple, Callable
from datetime import datetime
from calendar import monthrange
import requests
from requests.adapters import HTTPAdapter, Retry
import pandas as pd

logger = logging.getLogger(__name__)


def _get_month_range(year: int, month: int) -> Tuple[str, str]:
    """
    Returns the first and last day of a given month in a specified year.

    Args:
        year (int): The year.
        month (int): The month.

    Returns:
        Tuple[str, str]: A tuple containing the first and last day of the month.
    """
    first_day = datetime(year, month, 1)
    _, last_day_num = monthrange(year, month)
    last_day = datetime(year, month, last_day_num)
    return first_day.strftime("%Y-%m-%d"), last_day.strftime("%Y-%m-%d")


def _setup_session() -> requests.Session:
    """Set up the requests session with retry strategy."""
    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=0.3,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["POST"],
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session


def _get_nih_projects(
    session: requests.Session, api_url: str, query_params: Dict[str, str]
) -> Dict[str, str]:
    """
    Make an API POST request to fetch NIH projects.

    Args:
        api_url (str): API URL to use to fetch NIH funded projects.
        session (requests.Session): The requests session.
        query_params (dict): The query parameters for the API call.

    Returns:
        dict: The response data as a dictionary.
    """
    try:
        response = session.post(api_url, json=query_params, timeout=60)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        logger.error("Request failed: %s", e)


def _create_query_params(
    from_date: str, to_date: str, offset: int, limit: int
) -> Dict[str, str]:
    """
    Creates query parameters for NIH API requests.

    Args:
        from_date (str): Start date for the query range.
        to_date (str): End date for the query range.
        offset (int): The offset for pagination.
        limit (int): Maximum number of records per API request.

    Returns:
        dict: A dictionary of query parameters.
    """
    return {
        "criteria": {"award_notice_date": {"from_date": from_date, "to_date": to_date}},
        "offset": offset,
        "limit": limit,
    }


def _collect_nih_funded_projects_per_month(
    api_url: str,
    year: int,
    month: int,
    limit: int = 500,
) -> pd.DataFrame:
    """
    Fetches NIH funded projects for a specified year and month.

    Args:
        api_url (str): API URL to use to fetch NIH funded projects.
        year (int): The fiscal year for which to fetch the grants.
        month (int): The month for which to fetch the grants.
        limit (int): Maximum number of records per API request.

    Returns:
        pd.DataFrame: A DataFrame of funded projects for the specified year and month.
    """
    session = _setup_session()
    from_date, to_date = _get_month_range(year, month)

    all_projects = []
    offset = 0

    while True:
        time.sleep(1)
        query_params = _create_query_params(from_date, to_date, offset, limit)
        data = _get_nih_projects(session, api_url, query_params)
        if offset == 0:
            total_count = data["meta"]["total"]
        projects = data.get("results", [])
        all_projects.extend(projects)
        offset += len(projects)
        if offset % 1000 == 0 or len(projects) < limit:
            logger.info(
                "Collected %d/%d NIH funded projects for %d-%02d",
                offset,
                total_count,
                year,
                month,
            )
        if len(projects) < limit:
            break

    return pd.json_normalize(all_projects)


def _collect_nih_funded_projects_per_year(year: int, api_url: str) -> pd.DataFrame:
    """
    Fetches NIH funded projects for a specified year, month by month.

    Args:
        year (int): The year for which to fetch the grants.
        api_url (str): API URL to use to fetch NIH funded projects.

    Returns:
        pd.DataFrame: A DataFrame containing all funded projects for the specified year.
    """
    all_projects = []
    for month in range(1, 13):
        monthly_projects = _collect_nih_funded_projects_per_month(api_url, year, month)
        all_projects.append(monthly_projects)
    return pd.concat(all_projects, ignore_index=True)


def collect_nih_funded_projects(years: List[int], api_url: str) -> Dict[str, Callable]:
    """
    Generates a dictionary with keys as dataset names and values as callable
    functions to fetch NIH funded projects for each year in the provided list.

    Args:
        years (List[int]): A list of years for which to generate the datasets.
        api_url (str): API URL to use to fetch NIH funded projects.

    Returns:
        Dict[str, Callable]: A dictionary where each key is a dataset name for
            a specific year, and each value is a callable that fetches the dataset
            for that year.
    """
    return {
        f"nih_funded_projects_{year}": lambda year=year: _collect_nih_funded_projects_per_year(
            year, api_url
        )
        for year in years
    }
