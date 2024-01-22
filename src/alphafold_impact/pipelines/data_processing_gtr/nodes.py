"""
This module contains the functions used in the data processing pipeline for
    the Gateway to Research data.

The functions are used to load the data, process it, and save it in the
    required format.
"""
import logging
from typing import Dict, List, Tuple
import pandas as pd
from kedro.io import AbstractDataset
from joblib import Parallel, delayed

logger = logging.getLogger(__name__)


def load_institutions_data(data_dict: Dict[str, AbstractDataset]) -> Dict[str, dict]:
    """
    Load institutions data from a dictionary of partitions.

    Args:
        data_dict (dict): A dictionary containing partitions of data.

    Returns:
        dict: A dictionary containing the loaded institutions data.
    """
    logger.info("Loading chunked data.")
    gtr_data = Parallel(n_jobs=-1, verbose=10)(
        delayed(_load_institutions)(partition) for partition in data_dict.values()
    )

    # split the list of dictionaries into two lists of dictionaries
    gtr_authors, gtr_publications = zip(*gtr_data)

    # join the list of dictionaries into a single dictionary
    return (
        {k: v for d in gtr_authors for k, v in d.items()},
        {k: v for d in gtr_publications for k, v in d.items()},
    )


def _load_institutions(loader_obj: AbstractDataset) -> Dict[str, Tuple[str, str, str]]:
    """
    Load institutions data from the given loader object.

    Args:
        loader_obj (AbstractDataset): The loader object to load data from.

    Returns:
        dict: A dictionary containing the institutions data.
            The keys are work IDs, and the values are dictionaries
            where the keys are author IDs and the values are lists of
            institution information.
    """
    data = loader_obj()
    institutions_dict, publications_dict = {}, {}
    for chunk in data:
        for work in chunk:
            institutions_dict[work_id := work["id"].split("/")[-1]] = {}
            publications_dict[work_id] = {
                "title": work["title"],
                "abstract": work["abstract"],
                "doi": work["doi"],
                "publication_date": work["publication_date"],
            }
            for author in work["authorships"]:
                try:
                    institution_list = []
                    for institution in author.get("institutions", []):
                        triple = [
                            institution.get("id", "").split("/")[-1],
                            institution.get("display_name", ""),
                            institution.get("country_code", None),
                        ]
                        institution_list.append(triple)
                    institutions_dict[work_id][
                        author["author"]["id"].split("/")[-1]
                    ] = institution_list
                except KeyError:
                    pass
    return institutions_dict, publications_dict


def process_institutions(gtr_authors: List[dict]):
    """
    Process the institutions data from GTR authors.

    Args:
        gtr_authors (List[dict]): List of dictionaries containing GTR authors data.

    Returns:
        pd.DataFrame: Processed institutions data as a pandas DataFrame.
    """

    logger.info("Processing institutions data.")

    # unpack dictionaries, remove url from oa id
    gtr_authors = [
        (k, author_id, institution_id, institution_name, country_code)
        for k, v in gtr_authors.items()
        for author_id, institution_list in v.items()
        for institution_id, institution_name, country_code in institution_list
    ]

    # create a pandas dataframe
    data_institutions = pd.DataFrame(
        gtr_authors,
        columns=[
            "work_id",
            "author_id",
            "institution_id",
            "institution_name",
            "country_code",
        ],
    )

    return data_institutions

def process_publications(gtr_publications: List[dict]):
    """
    Process the publications data from GTR publications.

    Args:
        gtr_publications (List[dict]): List of dictionaries containing GTR publications data.

    Returns:
        pd.DataFrame: Processed publications data as a pandas DataFrame.
    """

    logger.info("Processing publications data.")

    # unpack dictionaries, remove url from oa id
    gtr_publications = [
        (k, v["title"], v["abstract"], v["doi"], v["publication_date"])
        for k, v in gtr_publications.items()
    ]

    # create a pandas dataframe
    data_publications = pd.DataFrame(
        gtr_publications,
        columns=[
            "work_id",
            "title",
            "abstract",
            "doi",
            "publication_date",
        ],
    )

    return data_publications


def process_alphafold_citations(reference_dict, cites_data):
    """
    Process the AlphaFold citations for each work in the reference dictionary.

    Args:
        reference_dict (dict): A dictionary containing works as keys and their
            corresponding information as values.
        cites_data (dict): A dictionary containing citation data.

    Returns:
        pandas.DataFrame: A dataframe containing the processed AlphaFold
            citations for each work.
    """

    # load alphafold
    af_data = cites_data["W3177828909"]()

    # unpack dictionaries, remove url from oa id
    list_af_cites = [v["id"].split("/")[-1] for v in af_data]

    # add af_direct_citation and af_indirect_citations to reference_dict
    for k, v in reference_dict.items():
        reference_dict[k]["af_indirect_citations"] = len(
            set(list_af_cites).intersection(v["referenced_works"])
        )

        reference_dict[k]["af_direct_citation"] = len(
            set(["W3177828909"]).intersection(v["referenced_works"])
        )

    # create a pandas dataframe for work_id, doi, af_direct_citation, af_indirect_citations
    af_citations = pd.DataFrame.from_dict(reference_dict, orient="index")
    af_citations = af_citations.reset_index()
    af_citations = af_citations.rename(columns={"index": "work_id"})
    af_citations = af_citations[
        ["work_id", "doi", "af_direct_citation", "af_indirect_citations"]
    ]
    return af_citations
