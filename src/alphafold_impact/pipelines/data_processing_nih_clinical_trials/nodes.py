"""
This module contains the functions used in the data processing pipeline for
the NIH clinical trials data.

Functions:
    trials_with_references: Loads NIH clinical trials from partitioned
        dataset, filters out trials that do not reference a paper,
        concatenates results into one DataFrame.
    nih_clinical_trials_links_to_papers: Produces a DataFrame which can
        be used to link NIH clinical trials to research papers via 
        PubMed ID or DOI.
"""
from typing import Dict, Callable, Any
import logging
import ast
import re
import pandas as pd


logger = logging.getLogger(__name__)


def trials_with_references(
    nih_clinical_trials_partitions: Dict[str, Callable[[], Any]]
) -> pd.DataFrame:
    """Loads NIH clinical trials from partitioned dataset, filters out trials
    that do not reference a paper, concatenates results into one DataFrame.

    Args:
        nih_clinical_trials_partitions (dict): Partitioned dataset of NIH clinical trials.

    Returns:
        DataFrame of NIH clinical trials with a paper reference.
    """
    ct_with_refs_combined = pd.DataFrame()
    n_partitions = len(nih_clinical_trials_partitions)

    for idx, (_, ct_partition_load_func) in enumerate(
        sorted(nih_clinical_trials_partitions.items()), start=1
    ):
        ct_partition_data_with_refs = ct_partition_load_func().query(
            "`ProtocolSection.ReferencesModule.ReferenceList.Reference`.notna()"
        )

        ct_with_refs_combined = pd.concat(
            [ct_with_refs_combined, ct_partition_data_with_refs],
            ignore_index=True,
            join="outer",
        )

        logger.info(
            "Concatenated clinical trials with references from partition %s/%s",
            idx,
            n_partitions,
        )

    return ct_with_refs_combined


def _reference_string_to_df(row: pd.Series) -> pd.DataFrame:
    """
    Processes a row of the DataFrame, converting the reference string to a list of dictionaries.
    Returns a DataFrame where each dictionary becomes a row, and each key in the dictionaries becomes a column.

    Args:
        row (pd.Series): A row of the DataFrame.

    Returns:
        pd.DataFrame: A DataFrame where each row corresponds to a dictionary in the reference list.
    """
    return pd.DataFrame(
        ast.literal_eval(
            row["ProtocolSection.ReferencesModule.ReferenceList.Reference"]
        )
    ).assign(nct_id=row["ProtocolSection.IdentificationModule.NCTId"])


def _extract_doi(citation: str) -> str:
    """
    Extracts the DOI from a citation string, ensuring not to include a trailing period.

    Args:
        citation (str): The citation string from which to extract the DOI.

    Returns:
        str: The extracted DOI, or an empty string if no DOI is found.
    """
    # This pattern looks for 'doi: ' followed by a sequence of characters
    doi_pattern = r"doi: ([\S]+)"
    if not (match := re.search(doi_pattern, citation)):
        return ""
    doi = match[1]
    # Remove the trailing period if it's the last character
    if doi.endswith("."):
        doi = doi[:-1]
    return doi


def trials_links_to_papers(
    clinical_trials_with_references: pd.DataFrame,
) -> pd.DataFrame:
    """
    Processes clinical trial data with references to papers, to produce a DataFrame
    which can be used to link NIH clinical trials to research papers via PubMed ID
    or DOI.

    This function takes a DataFrame containing clinical trial data with references
    to papers and performs the following data processing steps:
    1. Creates new columns for each component of the paper reference.
    2. Creates new column for the DOI.
    3. Renames and reorders columns.

    Args:
        clinical_trials_with_references (pd.DataFrame): A DataFrame containing
            clinical trial data with references to papers.

    Returns:
        pd.DataFrame: A processed DataFrame containing the following columns:
            - 'nct_id': The NCT ID of the clinical trial.
            - 'ref_pmid': The PMID (PubMed ID) of the paper reference.
            - 'ref_doi': The DOI (Digital Object Identifier) of the paper reference.
            - 'ref_citation': The citation of the paper reference.
            - 'ref_type': The type of the paper reference.
            - 'ref_retraction_list': The retraction list associated with the paper reference.
    """
    return (
        clinical_trials_with_references.apply(_reference_string_to_df, axis=1)
        .pipe(lambda x: pd.concat(x.tolist(), ignore_index=True))
        .assign(ref_doi=lambda x: x["ReferenceCitation"].apply(_extract_doi))
        .rename(
            columns={
                "ReferencePMID": "ref_pmid",
                "ReferenceType": "ref_type",
                "ReferenceCitation": "ref_citation",
                "RetractionList": "ref_retraction_list",
            }
        )[
            [
                "nct_id",
                "ref_pmid",
                "ref_doi",
                "ref_citation",
                "ref_type",
                "ref_retraction_list",
            ]
        ]
    )
