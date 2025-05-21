import logging
from typing import Sequence, Dict, Generator
import requests
from requests.adapters import HTTPAdapter, Retry
import gzip
import json
import pandas as pd

logger = logging.getLogger(__name__)


def get_s2_presigned_urls(
    bulk_url: str,
    api_key: str,
) -> Sequence[str]:
    """ """

    headers = {"X-API-KEY": api_key}

    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=5,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))

    datasets_response = session.get(
        bulk_url,
        headers=headers,
    )

    datasets_response.raise_for_status()

    datasets_metadata = datasets_response.json()

    logger.info(
        "Obtained file links. %s",
        datasets_metadata["description"],
    )

    # create a dict of file_id -> presigned_url
    presigned_urls = [
        {"id": f"s{i}", "url": file_url}
        for i, file_url in enumerate(datasets_metadata["files"])
    ]

    return presigned_urls


def get_s2_bulk_data(
    presigned_url: Sequence[Dict[str, str]],
) -> Generator[Dict[str, pd.DataFrame], None, None]:
    """Fetch bulk data from S2 presigned URLs, extract JSON and convert to DataFrames."""

    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=5,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))

    for i, url in enumerate(presigned_url):
        logger.info(
            "Fetching %s. Number %s of %s", url["id"], i + 1, len(presigned_url)
        )
        response = session.get(url["url"])
        response.raise_for_status()

        json_data = json.loads(gzip.decompress(response.content).decode("utf-8"))
        df = pd.DataFrame(json_data)
        yield {url["id"]: df}
