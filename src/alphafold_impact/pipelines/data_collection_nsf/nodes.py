"""
This module contains the functions used to fetch data from the NSF API.

Functions:
    _create_publication_objects: Create publication objects from a list of paper data.
    fetch_award_data: Fetch award data for a given award ID.
    fetch_nsf_data: Fetch NSF data for a given award year.
"""

import logging
from typing import Dict, List
import requests
from requests.adapters import HTTPAdapter, Retry
from bs4 import BeautifulSoup as bs
from joblib import Parallel, delayed
import pandas as pd

logger = logging.getLogger(__name__)


def _create_publication_objects(publication_list: List[str]) -> Dict[str, str]:
    """
    Create publication objects from a list of paper data.

    Args:
        publication_list (List[str]): A list of paper data in the format
            "journal~year~volume~authors~doi~...".

    Returns:
        Dict[str, str]: A dictionary of publication objects, where the keys
            are the NSF PAR IDs and the values are dictionaries containing
            the publication information (journal, year, volume, authors,
            title, doi).
    """
    publication_dict = {}
    for paper in publication_list:
        paper_data = paper.split("~")
        if len(paper_data) != 17:
            logger.info("Paper data is incomplete. Skipping paper.")
            continue
        else:
            journal = paper_data[0]
            year = paper_data[1]
            volume = paper_data[2]
            authors = paper_data[3]
            doi = paper_data[4]
            title = paper_data[9]
            nsf_par_id = paper_data[13]

            publication_dict[nsf_par_id] = {
                "journal": journal,
                "year": year,
                "volume": volume,
                "authors": authors,
                "title": title,
                "doi": doi,
            }

    return publication_dict


def fetch_award_data(
    config: Dict[str, str], fields: List[str], award_id: str
) -> Dict[str, str]:
    """
    Fetches award data for a given award ID.

    Args:
        config (Dict[str, str]): Configuration parameters.
        fields (List[str]): List of fields to include in the response.
        award_id (str): ID of the award to fetch data for.

    Returns:
        Dict[str, str]: Award data as a dictionary.

    """
    base_url = config["base_url"]
    fields = ",".join(fields)
    session = requests.Session()
    retries = Retry(
        total=config["max_retries"], backoff_factor=config["backoff_factor"]
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))

    logger.info("Requesting award data for %s", award_id)
    try:
        while True:
            response = session.get(
                f"{base_url}?id={award_id}&printFields={fields}",
                timeout=30,
            )
            if not response.content:
                logger.info("Received empty response for %s", award_id)

            if response.status_code == 200 and response.content:
                break

        data = response.json()
        award = data.get("response", {}).get("award", [])
        if award:
            logger.info("Received award data for %s", award_id)
            award_object = award[0]
            if "publicationResearch" in award_object:
                award_object["publicationResearch"] = _create_publication_objects(
                    award_object["publicationResearch"]
                )
            else:
                award_object["publicationResearch"] = {}
            return award_object
        else:
            logger.info("No award data received for %s", award_id)
            return {}
    except Exception as e:  # pylint: disable=broad-except
        logger.exception("Failed to fetch award data for %s: %s", award_id, e)
        return {}


def fetch_nsf_data(
    config: Dict[str, str], fields: List[str], award_list: List[str]
) -> Dict[str, Dict[str, str]]:
    """
    Fetches NSF data for a given list of awards.

    Args:
        config (Dict[str, str]): Configuration parameters for fetching data.
        fields (List[str]): List of fields to fetch for each award.
        award_list (str): List of award IDs to fetch data for.

    Returns:
        dict: A dictionary containing the fetched NSF data.
    """
    logger.info("Fetching NSF data for year %s", award_list[0][:2])
    awards = Parallel(n_jobs=5, verbose=10)(
        delayed(fetch_award_data)(config, fields, award_id) for award_id in award_list
    )
    logger.info("Fetched NSF data for year %s", award_list[0][:2])
    awards = [award for award in awards if award.get("id")]
    return {d["id"]: {k: v for k, v in d.items() if k != "id"} for d in awards}


def fetch_nsf_archived_funding_opportunities(
    archived_search_url: str,
    config: Dict[str, str],
    base_url: str,
) -> pd.DataFrame:
    """
    Fetches archived funding opportunities from the NSF website and returns 
        them as a pandas DataFrame.

    Args:
        archived_search_url (str): The URL of the archived search page.
        config (Dict[str, str]): Configuration parameters for the request.
        base_url (str): The base URL of the NSF website.

    Returns:
        pd.DataFrame: A DataFrame containing the fetched funding opportunities 
            with the following columns:
            - URL: The URL of the funding opportunity.
            - Title: The title of the funding opportunity.
            - Synopsis: The synopsis of the funding opportunity.
            - Linked Awards: The linked awards of the funding opportunity.
            - Solicitation URL: The URL of the solicitation.

    """
    session = requests.Session()
    retries = Retry(
        total=config["max_retries"], backoff_factor=config["backoff_factor"]
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))
    response = session.get(archived_search_url)

    soup = bs(response.content, "html.parser")
    opportunity_links = soup.find_all("p", class_="l-exc__heading")
    funding_opps = [
        {
            "href": base_url + link.a["href"],
            "title": link.a.get_text(strip=True).replace("\xa0", " "),
        }
        for link in opportunity_links
        if link.a
    ]

    # get opportunity details for each funding opportunity
    details = [
        _get_opportunity_details(funding_opp["href"], base_url, session)
        for funding_opp in funding_opps
    ]

    # add the details to the funding opportunities
    for i, funding_opp in enumerate(funding_opps):
        funding_opp.update(details[i])

    # transform to dataframe with custom column names
    return pd.DataFrame(
        funding_opps,
        columns=["URL", "Title", "Synopsis", "Linked Awards", "Solicitation URL"],
    )


def _get_opportunity_details(
    href: str, base_url: str, session: requests.Session
) -> Dict[str, str]:
    """
    Retrieves the program details from a given URL.

    Args:
        href (str): The URL to retrieve the program details from.
        base_url (str): The base URL used to construct the complete links.
        session (requests.Session): The session object to send the HTTP GET request.

    Returns:
        dict: A dictionary containing the program details, including the synopsis,
            funded link, and solicitation link.
    """

    logger.info("Fetching details for %s", href)
    # Send an HTTP GET request to the URL
    response = session.get(href)

    # Parse the HTML content
    soup = bs(response.text, "html.parser")

    # Find the SYNOPSIS section and collect subsequent <p> elements
    synopsis_tag = soup.find("strong", class_="greybold", string="SYNOPSIS")
    synopsis_text = []

    if synopsis_tag:
        try:
            synopsis_text = "\n".join(
                [
                    sibling.text.strip()
                    for sibling in synopsis_tag.find_next()
                    if sibling.name == "p"
                ]
            )

        except TypeError:
            synopsis_text = ""

    # Find the right links
    try:
        funded_link_href = soup.find(
            "a",
            href=lambda href: href
            and (
                href.startswith("/awardsearch/advancedSearchResult")
                or href.startswith(base_url + "/awardsearch/advancedSearchResult")
            ),
        )["href"]

        if not funded_link_href.startswith(base_url):
            funded_link_href = base_url + funded_link_href

        funded_link_href = funded_link_href.strip()
    except TypeError:
        funded_link_href = ""
        
    try:
        solicitation_link_href = (
            base_url
            + soup.find(
                "a",
                href=lambda href: href
                and (
                    href.startswith("/publications/pub_summ")
                    or href.startswith(base_url + "/publications/pub_summ")
                ),
            )["href"]
        ).strip()
    except TypeError:
        solicitation_link_href = ""

    logger.info("Fetched details for %s", href)
    return {
        "synopsis": synopsis_text,
        "funded_link_href": funded_link_href,
        "solicitation_link_href": solicitation_link_href,
    }
