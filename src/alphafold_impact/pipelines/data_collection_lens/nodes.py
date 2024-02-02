"""
This module contains the functions used to fetch data from the Lens API.

Functions:
    create_request_form: Create the JSON request body and headers for the Lens API.
    fetch_lens_data: Fetch data from Lens API using scroll requests.
"""
import logging
import time
import yaml
from typing import Dict, Tuple
from calendar import monthrange
import requests
from requests.adapters import HTTPAdapter, Retry

logger = logging.getLogger(__name__)

def get_app_credentials():
    """Get the credentials for the Lens API."""
    with open("./conf/base/credentials.yml") as cred:
        token = yaml.safe_load(cred).get("lens_token")
    return token

def create_request_form(
    config: Dict[str, str], month: str, year: str, jurisdiction: str, token: str
) -> Tuple[str, Dict[str, str]]:
    """Create the JSON request body and headers for the Lens API.

    Args:
        config (Dict[str, str]): Configuration parameters.
        month (str): The month to query.
        year (str): The year to query.
        jurisdiction (str): The jurisdiction to query.

    Returns:
        Tuple[str, Dict[str, str]]: The JSON request body and headers.
    """
    start = f"{year}-{month}-01"
    end = f"{year}-{month}-{monthrange(int(year), int(month))[1]}"

    request_body = """{
        "query": {
            "bool": {
                "must": [
                    {
                        "match": {
                            "jurisdiction": "%s"
                        }
                    }, 
                    {
                        "range": {
                            "date_published": {
                                "gte": "%s", 
                                "lte": "%s"
                            }
                        }
                    }
                ]
            }
        },
        "size": %d,
        "scroll": "1m"
    }""" % (
        jurisdiction,
        start,
        end,
        config["size"],
    )

    headers = {"Authorization": token, "Content-Type": "application/json"}

    return request_body, headers


def fetch_lens_data(
    config: Dict[str, str], request_body: Dict[str, str], headers: Dict[str, str]
) -> Dict[str, str]:
    """Fetch data from Lens API using scroll requests.

    Args:
        config (Dict[str, str]): Configuration parameters.
        request_body (Dict[str, str]): The JSON request body.
        headers (Dict[str, str]): The request headers.

    Returns:
        Dict[str, str]: A dictionary of patent data. Keys are patent IDs,
            values are patent data.
    """
    lens_list = []
    session = requests.Session()
    retries = Retry(
        total=config["max_retries"], backoff_factor=config["backoff_factor"]
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))
    logger.info("Making first request to Lens API")
    base_url = config["base_url"]
    response = session.post(
        base_url,
        data=request_body,
        headers=headers,
        timeout=30,
    )

    if response.status_code == 200:
        first_response = response.json()
        lens_list.append(first_response["data"])
        logger.info("First response received.")
        logger.info("Total patents: %s", first_response["total"])

    if "scroll_id" in first_response:
        scroll_request = {"scroll_id": first_response["scroll_id"], "scroll": "1m"}

        while True:
            time.sleep(8)  # per minute rate limit == 9
            logger.info("Making scroll request to Lens API")
            response = session.post(
                base_url,
                json=scroll_request,
                headers=headers,
                timeout=30,
            )
            if response.status_code == 200:
                scroll_response = response.json()
                lens_list.append(scroll_response["data"])
                logger.info("Scroll response received.")
            elif response.status_code == 400:
                logger.info("Scroll request exhausted, no response.")
                break
            elif response.status_code == 429:
                logger.info("Scroll request rate limit exceeded.")
                time.sleep(60)
                continue
            else:
                logger.info(
                    "Scroll request failed. Response code: %s", response.status_code
                )
                break

    # Flatten list of lists
    lens_list = [item for sublist in lens_list for item in sublist]

    # Refactor as a dictionary, with "lens_id" as the key
    lens_dict = {item["lens_id"]: item for item in lens_list}

    return lens_dict
