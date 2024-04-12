"""Pipeline nodes for the data_collection_s2 pipeline.

Functions:
    fetch_citation_details: Fetches citation details for a given work ID.
    process_citations: Process citation outputs and extract relevant 
        information.
    process_references: Process the reference outputs and generate a 
        list of rows containing relevant information.
    iterate_citation_detail_points: Iterates over citation detail 
        points and fetches citation details based on the provided parameters.
    get_intent_level_0: Retrieves the intent level 0 data from the given 
        OA dataset.
    get_intent_level: Retrieves the intent level data from the given 
        OA dataset based on the specified level.
    get_citation_intent_from_oa_dataset: Retrieves citation intent data 
        from an Open Access dataset.
"""

import logging
from typing import Sequence, Dict, Union, List, Any
import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from joblib import Parallel, delayed

logger = logging.getLogger(__name__)


def fetch_citation_details(
    work_id: str,
    base_url: str,
    direction: str,
    fields: Sequence[str],
    api_key: str,
    perpage: int = 500,
) -> Sequence[Dict[str, str]]:
    """
    Fetches citation details for a given work ID.

    Args:
        work_id (str): The work ID to fetch citation details for.
        base_url (str): The base URL for the API.
        direction (str): The direction of the citations.
        fields (Sequence[str]): The fields to fetch.
        perpage (int, optional): The number of citations to fetch per page.
            Defaults to 500.

    Returns:
        Sequence[Dict[str, str]]: A list of citation details.
    """
    offset = 0
    url = (
        f"{base_url}/{work_id}/{direction}?"
        f"fields={','.join(fields)}&offset={offset}&limit={perpage}"
    )

    headers = {"X-API-KEY": api_key}

    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=0.3,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))

    data_list = []

    while True:
        response = session.get(url, headers=headers)
        response.raise_for_status()
        data = response.json().get("data", [])
        data_list.extend(data)

        if len(data) < perpage:
            break

        offset += perpage
        url = (
            f"{base_url}/{work_id}/{direction}?"
            f"fields={','.join(fields)}&offset={offset}&limit={perpage}"
        )

    return data_list


# Updated functions
def process_citations(citation_outputs: Dict[str, Sequence[Dict[str, Any]]]) -> List:
    """
    Process citation outputs and extract relevant information.

    Args:
        citation_outputs (Dict[str, Sequence[Dict[str, Any]]]): A dictionary
            containing citation outputs.

    Returns:
        List: A list of rows containing extracted information from the
            citation outputs.
    """
    rows = []
    for parent_doi, outputs in citation_outputs.items():
        for output in outputs:
            context_intents = output.get("contextsWithIntent", [])
            influential = output.get("isInfluential", "")
            external_ids = output.get("citingPaper", {}).get("externalIds", {})
            if not external_ids:
                continue
            doi = external_ids.get("DOI", "")
            pmid = external_ids.get("PubMed", "")
            for item in context_intents:
                context = item.get("context", "")
                intents = item.get("intents") or [""]
                for intent in intents:
                    rows.append(
                        [
                            parent_doi,
                            pmid,
                            doi,
                            influential,
                            intent,
                            context,
                        ]
                    )

    return rows


def process_references(reference_outputs: Dict[str, Dict[str, Any]]) -> List:
    """
    Process the reference outputs and generate a list of rows containing relevant
        information.

    Args:
        reference_outputs (Dict[str, Dict[str, Any]]): A dictionary containing the
            reference outputs.

    Returns:
        List: A list of rows containing the processed information.

    """
    rows = []
    for child_doi, output in reference_outputs.items():
        context_intents = output.get("contextsWithIntent", [])
        influential = output.get("isInfluential", "")
        parent_doi = output.get("citedPaper", {}).get("externalIds", {}).get("DOI", "")
        parent_pmid = (
            output.get("citedPaper", {}).get("externalIds", {}).get("PubMed", "")
        )
        for item in context_intents:
            context = item.get("context", "")
            intents = item.get("intents") or [""]
            for intent in intents:
                rows.append(
                    [
                        parent_doi,
                        parent_pmid,
                        child_doi,
                        influential,
                        intent,
                        context,
                    ]
                )
    return rows


def iterate_citation_detail_points(
    oa: str,
    parent_doi: str,
    doi: str,
    parent_pmid: str,
    pmid: str,
    direction: str,
    **kwargs,
) -> Dict[str, Union[str, Sequence[Dict[str, Any]]]]:
    """
    Iterates over citation detail points and fetches citation details based on the
        provided parameters.

    Args:
        oa (str): The OA (Open Access) identifier.
        parent_doi (str): The DOI (Digital Object Identifier) of the parent paper.
        doi (str): The DOI of the current paper.
        parent_pmid (str): The PubMed ID of the parent paper.
        pmid (str): The PubMed ID of the current paper.
        direction (str): The direction of the citation details to fetch. Can be
            "references" or "forward".
        **kwargs: Additional keyword arguments to pass to the fetch_citation_details
            function.

    Returns:
        Dict[str, Union[str, Sequence[Dict[str, Any]]]]: A dictionary containing the fetched
            citation details.

    """
    logger.info("Fetching citation details for %s", oa)
    for prefix, id_ in [("DOI:", doi), ("PMID:", pmid)]:
        if not id_:
            logger.warning("No relevant %s id", prefix[:-1])
            continue
        work_id = f"{prefix}{id_}"
        try:
            data = fetch_citation_details(
                work_id=work_id, direction=direction, **kwargs
            )

            # if looking for a relevant backwards citation, subset to parent_oa
            if direction == "references":
                for item in data:
                    if (eids := item["citedPaper"]["externalIds"]) is not None:
                        if (
                            eids.get("DOI", "") == parent_doi
                            or eids.get("PubMed", "") == parent_pmid
                        ):
                            logger.info("Found citation details using %s", work_id)
                            return item
                logger.warning("No relevant citation found in %s", work_id)

            # otherwise, direction is all forward citations, keep all details
            else:
                logger.info("Found citation details using %s", work_id)
                return data
        except requests.exceptions.HTTPError as errh:
            logger.error("HTTP Error: %s", errh)
            return {}
        except requests.exceptions.ConnectionError as errc:
            logger.error("Error Connecting: %s", errc)
            return {}
        except requests.exceptions.Timeout as errt:
            logger.error("Timeout Error: %s", errt)
            return {}
        except requests.exceptions.RequestException as err:
            logger.error("Something went wrong: %s", err)
            return {}
        except Exception:  # pylint: disable=broad-except
            logger.warning("Failed to fetch citation details for %s", work_id)
            return {}
    return {}


def get_intent_level_0(oa_dataset: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """
    Retrieves the intent level 0 data from the given OA dataset. Note that the
    pipeline is different for level 0. externalIds correspond to the upstream citation
    (reference), while in subsequent levels it corresponds to downstream citations.

    Args:
        oa_dataset (pd.DataFrame): The input OA dataset.
        **kwargs: Additional keyword arguments.

    Returns:
        pd.DataFrame: The processed intent level 0 data.
            - parent_doi: The DOI of the parent publication.
            - parent_pmid: The PMID of the parent publication.
            - doi: The DOI of the citation.
            - influential: Indicates whether the citation is influential.
            - intent: The intent of the citation.
            - context: The context of the citation.

    """
    # Let's do level zero separately as the pipeline is different
    inputs = (
        oa_dataset[oa_dataset["level"] == 0]
        .apply(
            lambda x: (x["id"], x["parent_doi"], x["doi"], x["parent_pmid"], x["pmid"]),
            axis=1,
        )
        .tolist()
    )

    level_outputs = Parallel(n_jobs=8)(
        delayed(iterate_citation_detail_points)(
            *input, direction="references", **kwargs
        )
        for input in inputs
    )

    # zip the child_oa with the outputs
    level_dict = dict(
        list(zip(oa_dataset[oa_dataset["level"] == 0]["doi"], level_outputs))
    )

    processed_references = process_references(level_dict)

    return pd.DataFrame(
        processed_references,
        columns=[
            "parent_doi",
            "parent_pmid",
            "doi",
            "influential",
            "intent",
            "context",
        ],
    )


def get_intent_level(oa_dataset: pd.DataFrame, level: int, **kwargs) -> pd.DataFrame:
    """
    Retrieves the intent level data from the given OA dataset based on the specified level.

    Args:
        oa_dataset (pd.DataFrame): The input OA dataset.
        level (int): The intent level to retrieve.
        **kwargs: Additional keyword arguments.

    Returns:
        pd.DataFrame: A DataFrame containing the processed intent level citations
            with the following columns:
            - parent_doi: The DOI of the parent publication.
            - pmid: The PMID of the citation.
            - doi: The DOI of the citation.
            - influential: Indicates whether the citation is influential.
            - intent: The intent of the citation.
            - context: The context of the citation.
    """

    level_data = oa_dataset[oa_dataset["level"] == level]

    # Remove duplicate rows based on "parent_id"
    level_data = level_data.drop_duplicates(subset="parent_id")

    inputs = level_data.apply(
        lambda x: (x["parent_id"], "", x["parent_doi"], "", x["parent_pmid"]),
        axis=1,
    ).tolist()

    # Use joblib to parallelize the function calls
    level_outputs = Parallel(n_jobs=8)(
        delayed(iterate_citation_detail_points)(*input, direction="citations", **kwargs)
        for input in inputs
    )

    level_dict = dict(list(zip(level_data["parent_doi"], level_outputs)))

    processed_level_citations = process_citations(level_dict)

    return pd.DataFrame(
        processed_level_citations,
        columns=[
            "parent_doi",
            "pmid",
            "doi",
            "influential",
            "intent",
            "context",
        ],
    )


def get_citation_intent_from_oa_dataset(
    oa_dataset: pd.DataFrame,
    **kwargs,
) -> pd.DataFrame:
    """
    Retrieves citation intent data from an Open Access dataset.

    Args:
        oa_dataset (pd.DataFrame): The Open Access dataset containing citation
            information.
        **kwargs: Additional keyword arguments.

    Returns:
        pd.DataFrame: The processed dataset with citation intent information.

    """
    # Create a mapping from id to doi and pmid for each level
    oa_dataset["parent_level"] = oa_dataset["level"] - 1

    oa_dataset = pd.merge(
        oa_dataset,
        oa_dataset,
        left_on=["parent_id", "parent_level"],
        right_on=["id", "level"],
        how="left",
        suffixes=("", "_parent"),
    )
    oa_dataset.rename(
        columns={"doi_parent": "parent_doi", "pmid_parent": "parent_pmid"}, inplace=True
    )

    # drop the other _parent columns
    oa_dataset.drop(
        columns=[col for col in oa_dataset.columns if "_parent" in col], inplace=True
    )

    # manually set the parent_pmid and parent_doi for AlphaFold
    oa_dataset.loc[oa_dataset["level"] == 0, "parent_doi"] = (
        "10.1038/s41586-021-03819-2"
    )
    oa_dataset.loc[oa_dataset["level"] == 0, "parent_pmid"] = "34265844"

    # get the citation intent for each level
    level_0 = get_intent_level_0(oa_dataset, **kwargs)
    level_1 = get_intent_level(oa_dataset, 1, **kwargs)
    level_2 = get_intent_level(oa_dataset, 2, **kwargs)
    level_3 = get_intent_level(oa_dataset, 3, **kwargs)

    # concatenate the dataframes
    processed_df = pd.concat([level_0, level_1, level_2, level_3])

    # check how often doi is empty and pmid is not
    logger.info(
        "Number of rows with empty doi and non-empty pmid: %d",
        processed_df[(processed_df["doi"] == "") & (processed_df["pmid"] != "")].shape[
            0
        ],
    )

    processed_grouped_data_doi = (
        processed_df.groupby(["parent_doi", "doi"])[
            ["influential", "intent", "context"]
        ]
        .apply(lambda x: x.to_dict(orient="records"))
        .reset_index()
        .rename(columns={0: "strength"})
    )

    processed_grouped_data_pmid = (
        processed_df.groupby(["parent_doi", "pmid"])[
            ["influential", "intent", "context"]
        ]
        .apply(lambda x: x.to_dict(orient="records"))
        .reset_index()
        .rename(columns={0: "strength"})
    )

    processed_data = pd.merge(
        oa_dataset, processed_grouped_data_doi, on=["parent_doi", "doi"], how="left"
    )

    processed_data = pd.merge(
        processed_data,
        processed_grouped_data_pmid,
        left_on=["parent_doi", "parent_pmid"],
        right_on=["parent_doi", "pmid"],
        how="left",
    )

    # combine the two strength columns
    processed_data["strength"] = processed_data["strength_x"].combine_first(
        processed_data["strength_y"]
    )

    # rename the pmid columns
    processed_data.rename(columns={"pmid_x": "pmid"}, inplace=True)
    processed_data.drop(columns=["strength_x", "strength_y", "pmid_y"], inplace=True)

    return processed_data[
        [
            "parent_id",
            "parent_doi",
            "parent_pmid",
            "id",
            "doi",
            "pmid",
            "level",
            "publication_date",
            "mesh_terms",
            "cited_by_count",
            "authorships",
            "parent_level",
            "strength",
        ]
    ]
