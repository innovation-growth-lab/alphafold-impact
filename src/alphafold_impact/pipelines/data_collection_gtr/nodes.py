"""
This module contains functions for fetching and preprocessing data from the GtR API.

Functions:
- fetch_gtr_data(parameters: Dict[str, Union[str, int]]) -> List[Dict[str, Any]]:
    Fetches data from the GtR API based on the provided parameters.

- preprocess_data_to_df(raw_data: List[Dict[str, Any]]) -> pd.DataFrame:
    Preprocesses the raw data into a pandas DataFrame.

- preprocess_organisations(org_df: pd.DataFrame) -> pd.DataFrame:
    Preprocesses the data by extracting the main address and dropping unnecessary columns.
"""

import logging
import random
import time
from typing import Dict, List, Union, Any
import requests
from requests.adapters import HTTPAdapter, Retry
import pandas as pd

logger = logging.getLogger(__name__)


def _extract_main_address(addresses):
    address = addresses["address"]
    main_address = next(
        (addr for addr in address if addr["type"] == "MAIN_ADDRESS"), None
    )
    return (
        pd.Series(main_address)
        if main_address
        else pd.Series([None] * 4, index=["id", "created", "postCode", "region"])
    )


def fetch_gtr_data(parameters: Dict[str, Union[str, int]]) -> List[Dict[str, Any]]:
    """Fetch data from the GtR API.

    Args:
        parameters (Dict[str, Union[str, int]]): Parameters for the API request.

    Returns:
        List[Dict[str, Any]]: The fetched data.
    """
    gtr_config = parameters["gtr_config"]
    endpoint, key = ( # [TEMP] Hardcoded for now on orgs
        gtr_config["orgs"]["endpoint"], 
        gtr_config["orgs"]["key"]
    )
    base_url, headers, page_size = (
        gtr_config["base_url"],
        gtr_config["headers"],
        gtr_config["page_size"]
    )
    max_retries, backoff_factor = (
        parameters["max_retries"],
        parameters["backoff_factor"]
    )

    all_data = []
    page = 1

    while True:
        url = f"{base_url}{endpoint}?page={page}&size={page_size}"
        s = requests.Session()
        retries = Retry(total=max_retries, backoff_factor=backoff_factor)
        s.mount("https://", HTTPAdapter(max_retries=retries))
        response = s.get(url, headers=headers, timeout=30)
        data = response.json()
        if key in data:
            items = data[key]
            if not items:
                break
            for item in items:
                item["page_fetched_from"] = page  # Add page info
                all_data.append(item)
        else:
            logger.error("No '%s' key found in the response", key)
            break

        logger.info("Fetched page %s from %s", page, endpoint)
        page += 1
        time.sleep(random.randint(3,10))  # [HACK] Respect web etiquette

    return all_data


def preprocess_data_to_df(raw_data: List[Dict[str, Any]]) -> pd.DataFrame:
    """Preprocess data to a DataFrame.

    Args:
        raw_data (List[Dict[str, Any]]): The raw data in a list of dictionaries.

    Returns:
        pd.DataFrame: The preprocessed data.
    """
    return pd.DataFrame(raw_data)


def preprocess_organisations(org_df: pd.DataFrame) -> pd.DataFrame:
    """Preprocess the organisations data.

    It extracts the main address and drops the "links" column.

    Args:
        org_df (pd.DataFrame): The organisations data.

    Returns:
        pd.DataFrame: The preprocessed data.
    """
    address_columns = org_df["addresses"].apply(_extract_main_address)
    address_columns = address_columns.drop("created", axis=1).add_prefix("address_")
    org_df = org_df.drop("addresses", axis=1).join(address_columns)
    org_df = org_df.drop(columns=["links"])
    return org_df
