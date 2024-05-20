"""
This module contains functions for
calculaing time to clinical citation (TCC)
"""

import logging
from urllib.error import HTTPError
import time
from typing import List, Union
import requests
import numpy as np
from Bio import Entrez
import pandas as pd

logger = logging.getLogger(__name__)


def convert_single_pubmed_id_to_int_string(
    value: Union[str, float]
) -> Union[str, float]:
    """
    Converts a value in the 'cited_by_clin' column to the correct format.

    If the value is a single number (float), it's converted to an int then to a string.
    If the value is a string with multiple numbers, it's left as is.
    NaN values are left unchanged.

    Args:
        value (float, str, or nan): The original value in the 'cited_by_clin' column.

    Returns:
        str or nan: The converted value.
    """
    # Check if the value is a single float (a single PubMed ID)
    if isinstance(value, float) and not pd.isna(value):
        return str(int(value))
    # Return the value unchanged if it's a string (multiple PubMed IDs) or NaN
    return value


def efetch_pmid(db: str, id: str, retmode: str, max_retries: int = 3, delay: int = 5):
    """
    Performs a Entrez efetch request to get DOI using PMID.

    Args:
        db (str): Database to query.
        id (str): The identifier of the article.
        retmode (str): Return mode for the request.
        max_retries (int): Maximum number of retries.
        delay (int): Delay between retries in seconds.

    Returns:
        HTTP response object on success, None on failure.
    """
    for attempt in range(max_retries):
        try:
            return Entrez.efetch(db=db, id=id, retmode=retmode)
        except HTTPError as e:
            if e.code == 502:
                logger.info(
                    "Attempt %d/%d for PMID %s failed with 502 Bad Gateway. Retrying in %d seconds...",
                    attempt + 1,
                    max_retries,
                    id,
                    delay,
                )
                time.sleep(delay)
            else:
                logger.error("HTTPError %d for PMID %s: %s", e.code, id, e.reason)
                return None
        except Exception as e:
            logger.error("Error fetching data for PMID %s: %s", id, str(e))
            return None

    logger.info("Maximum retries reached. Unable to fetch data for PMID %s.", id)
    return None


def get_doi_from_pubmed(pmid: str) -> str:
    """
    Fetches the DOI of a PubMed article using its PMID. If unable to fetch data,
    prints the PMID and the reason.

    Args:
        pmid (str): The PubMed ID of the article.

    Returns:
        str: The DOI of the article, or np.nan if DOI not found or an error occurs.
    """
    handle = efetch_pmid(db="pubmed", id=pmid, retmode="xml")
    if handle is None:
        return np.nan

    try:
        records = Entrez.read(handle)
    except Exception as e:
        logger.error("Error reading data for PMID %s: %s", pmid, str(e))
        return np.nan
    finally:
        handle.close()

    article = records["PubmedArticle"][0]["MedlineCitation"]["Article"]
    if "ELocationID" in article:
        for elocation in article["ELocationID"]:
            if elocation.attributes.get("EIdType") == "doi":
                return str(elocation)

    logger.info("DOI not found for PMID %s.", pmid)
    return np.nan


def fetch_pub_date_from_openalex(doi: str) -> str:
    """
    Fetches the publication date from OpenAlex using the DOI of the article.

    Args:
        doi (str): The DOI of the article.

    Returns:
        str: The publication date in YYYY-MM-DD format, or np.nan.
    """
    try:
        response = requests.get(f"https://api.openalex.org/works/doi:{doi}")
        data = response.json()
        publication_date = data.get("publication_date", None)
        if publication_date:
            return publication_date
        else:
            return np.nan
    except Exception:
        return np.nan


def get_publication_date(pmid: str) -> str:
    """
    Fetches the publication date of an article using its PMID.
    Tries PubMed first, then uses the DOI to fetch from OpenAlex.

    Args:
        pmid (str): The PubMed ID of the article.

    Returns:
        str: The publication date in YYYY-MM-DD format.
    """
    doi = get_doi_from_pubmed(pmid)
    return fetch_pub_date_from_openalex(doi)


def get_pub_dates_for_pmids(pmids: str) -> List[str]:
    """
    Fetches publication dates for a list of PubMed IDs.
    Skips if the input is NaN or non-string.

    Args:
        pmids (str): A string containing PubMed IDs separated by spaces.

    Returns:
        List[str]: A list of publication dates in YYYY-MM-DD format or np.nan.
                  Returns an empty list if pmids is NaN or non-string.
    """
    if pd.isna(pmids) or not isinstance(pmids, str):
        return np.nan
    return [get_publication_date(pmid) for pmid in pmids.split()]


def find_earliest_date(dates: Union[List[str], float]) -> Union[str, float]:
    """
    Finds the earliest date in a list of dates.

    Args:
        dates (Union[List[str], float]): A list of date strings or nan.

    Returns:
        Union[str, float]: The earliest date in string format or nan.
    """
    if isinstance(dates, list):
        dates = [date for date in dates if pd.notna(date)]
    else:
        return np.nan
    if dates == []:
        return np.nan
    datetime_list = pd.to_datetime(dates, errors="coerce")
    min_date = min(datetime_list)
    return min_date.strftime("%Y-%m-%d") if pd.notna(min_date) else np.nan


def match_oa_papers_to_icite(
    icite_data: dict,
    oa_data: pd.DataFrame,
    dataset_name: str,
) -> pd.DataFrame:
    """
    Find papers that are in iCite data based on both DOI and PMID.

    Args:
        icite_data (pd.DataFrame): iCite data with 'doi' and 'pmid' columns.
        oa_data (pd.DataFrame): OpenAlex data with 'doi' and 'pmid' columns.
        dataset_name (str): A label used for logging purposes, to denote the dataset
            being processed.

    Returns:
        pd.DataFrame: A DataFrame containing merged information from iCite and OpenAlex datasets.
                      It includes all columns from iCite and OpenAlex datasets.
    """
    oa_dois = oa_data["doi"].dropna().to_list()
    oa_pmids = [int(pmid) for pmid in oa_data.pmid.dropna().to_list()]

    logger.info("%s -- Finding matches in iCite data", dataset_name)

    doi_matches = icite_data.query(f"doi in {oa_dois}")
    pmid_matches = icite_data.query(f"pmid in {oa_pmids}")

    icite_matches = (
        pd.concat([doi_matches, pmid_matches])
        .drop_duplicates(subset="doi")
        .drop_duplicates(subset="pmid")
    )
    logger.info(
        "%s -- %0.2f%% OA papers matched to iCite (%d / %d)",
        dataset_name,
        len(icite_matches) / len(oa_data) * 100,
        len(icite_matches),
        len(oa_data),
    )
    return icite_matches


def add_clinal_citation_cols(icite_matches: pd.DataFrame) -> pd.DataFrame:
    """
    Modifies the provided DataFrame by converting 'cited_by_clin' values to integer strings
    and adding publication dates for the clinical articles.

    Args:
        icite_matches (pd.DataFrame): The DataFrame to be modified.

    Returns:
        pd.DataFrame: The modified DataFrame with updated 'cited_by_clin' and
                      'cited_by_clin_pub_dates' columns.
    """
    # Convert 'cited_by_clin' values to strings
    icite_matches["cited_by_clin"] = icite_matches["cited_by_clin"].apply(
        convert_single_pubmed_id_to_int_string
    )
    # Add number of clinical citations
    icite_matches["num_cited_by_clin"] = (
        icite_matches.cited_by_clin.str.count(" ").add(1).fillna(0, downcast="infer")
    )
    # Add clinical article publication dates
    icite_matches["cited_by_clin_pub_dates"] = icite_matches["cited_by_clin"].apply(
        get_pub_dates_for_pmids
    )
    # Add clinical article earliest publication date
    icite_matches["cited_by_clin_earliest_pub_date"] = icite_matches[
        "cited_by_clin_pub_dates"
    ].apply(find_earliest_date)
    icite_matches["cited_by_clin_earliest_pub_date"] = pd.to_datetime(
        icite_matches["cited_by_clin_earliest_pub_date"]
    )
    return icite_matches


def add_ctc_tcc(
    icite_data: pd.DataFrame,
    oa_data: pd.DataFrame,
    dataset_name: str,
    cutoff_date: str = "2021-07-15",
):
    """
    Combines OpenAlex data with iCite data based on DOI and PMID identifiers.
    Adds several columns related to clinical citations to the dataset,
    such as the number of clinical citations and the days until a clinical citation.

    Args:
        icite_data (pd.DataFrame): A dictionary containing iCite data with DOI and
            PMID keys.
        oa_data (pd.DataFrame): A DataFrame containing OpenAlex data with
            DOI and PMID columns.
        dataset_name (str): A string identifier used to name the dataset.
        cutoff_date (str): Tthe cutoff date for filtering publications.
            Defaults to "2021-07-15".

    Returns:
        pd.DataFrame: A DataFrame containing merged iCite and OpenAlex data with additional columns:
                      - 'num_cited_by_clin' (times cited by clinical article)
                      - 'cited_by_clin_pub_dates' (publication dates for clinical citations)
                      - 'cited_by_clin_earliest_pub_date' (earliest clinical citation date)
                      - 'days_to_clinical_trial' (days from publication to earliest clinical citation)
    """
    oa_data_filtered = oa_data.query(f"publication_date >= '{cutoff_date}'")
    oa_records_filtered = len(oa_data) - len(oa_data_filtered)
    logger.info("%s -- OpenAlex papers: %d", dataset_name, len(oa_data))
    logger.info(
        "%s -- OpenAlex papers after %s: %d",
        dataset_name,
        cutoff_date,
        len(oa_data_filtered),
    )
    logger.info(
        "%s -- OpenAlex papers removed after %s: %d",
        dataset_name,
        cutoff_date,
        oa_records_filtered,
    )
    icite_matches = match_oa_papers_to_icite(
        icite_data=icite_data,
        oa_data=oa_data_filtered,
        dataset_name=dataset_name,
    ).pipe(add_clinal_citation_cols)
    logger.info(
        "%s -- Number of clinical citations: %d",
        dataset_name,
        icite_matches.num_cited_by_clin.sum(),
    )
    logger.info(
        "%s -- Number of papers with clinical citations: %d",
        dataset_name,
        len(icite_matches.query("cited_by_clin.notna()")),
    )
    icite_oa_comb = icite_matches.merge(
        oa_data_filtered.drop_duplicates(subset="doi"),
        how="left",
        left_on="doi",
        right_on="doi",
    )
    icite_oa_comb["publication_date"] = pd.to_datetime(
        icite_oa_comb["publication_date"]
    )
    icite_oa_comb["days_to_clinical_trial"] = (
        icite_oa_comb["cited_by_clin_earliest_pub_date"]
        - icite_oa_comb["publication_date"]
    ).dt.days

    median_days_tcc = round(
        icite_oa_comb.query(
            "days_to_clinical_trial.notna()"
        ).days_to_clinical_trial.median(),
        1,
    )
    logger.info("%s -- Median days TCC: %f", dataset_name, median_days_tcc)
    return icite_oa_comb


def make_tcc_datasets(
    icite_data: pd.DataFrame,
    af_oa_data: pd.DataFrame,
    biochemistry_oa_data: pd.DataFrame,
    bioinformatics_oa_data: pd.DataFrame,
    protein_design_oa_data: pd.DataFrame,
    structural_biology_oa_data: pd.DataFrame,
    af_label: str,
    biochemistry_label: str,
    bioinformatics_label: str,
    protein_design_label: str,
    structural_biology_label: str,
):
    """
    For various subfields, combines data from iCite and OpenAlex based
    on DOI and PMID identifiers. Calculates time to clinical citation (TCC).

    Args:
        icite_data (pd.DataFrame): iCite data with 'doi' and 'pmid' columns.
        af_oa_data (pd.DataFrame): OpenAlex data for the papers that cite AF.
        biochemistry_oa_data (pd.DataFrame): OpenAlex data for the field of biochemistry.
        bioinformatics_oa_data (pd.DataFrame): OpenAlex data for the field of bioinformatics.
        protein_design_oa_data (pd.DataFrame): OpenAlex data for the field of protein design.
        structural_biology_oa_data (pd.DataFrame): OpenAlex data for the field of structural biology.
        af_label (str): Label for the papers that cite AF dataset.
        biochemistry_label (str): Label for the biochemistry dataset.
        bioinformatics_label (str): Label for the bioinformatics dataset.
        protein_design_label (str): Label for the protein design dataset.
        structural_biology_label (str): Label for the structural biology dataset.

    Returns:
        dict: A dictionary mapping each dataset label to a lambda function that, when called, will
              process and return a DataFrame containing the iCite and OpenAlex data with TCC for
              the respective field.
    """
    return {
        af_label: lambda: add_ctc_tcc(
            icite_data=icite_data, oa_data=af_oa_data, dataset_name=af_label
        ),
        biochemistry_label: lambda: add_ctc_tcc(
            icite_data=icite_data,
            oa_data=biochemistry_oa_data,
            dataset_name=biochemistry_label,
        ),
        bioinformatics_label: lambda: add_ctc_tcc(
            icite_data=icite_data,
            oa_data=bioinformatics_oa_data,
            dataset_name=bioinformatics_label,
        ),
        protein_design_label: lambda: add_ctc_tcc(
            icite_data=icite_data,
            oa_data=protein_design_oa_data,
            dataset_name=protein_design_label,
        ),
        structural_biology_label: lambda: add_ctc_tcc(
            icite_data=icite_data,
            oa_data=structural_biology_oa_data,
            dataset_name=structural_biology_label,
        ),
    }
