"""This module contains functions for fetching citation data to a specific depth and enriching
the data with additional MeSH tags.

The enrichment process is done in the following steps:
    1. Load a batch of data from the input dataset.
    2. Load JSON data, transform it into a DataFrame, and perform data transformations.
    3. Fetch additional mesh terms for DOIs in a DataFrame using Biopython's Entrez module.
        - Filter DataFrame for rows where mesh_terms is empty.
        - Create list of DOIs.
        - Fetch PubMed IDs for each DOI.
        - Use PubMed IDs to fetch mesh terms from full-data objects.
        - Update the DataFrame with the new mesh terms for the corresponding DOI, 
            removing the [DOI] part.

OpenAlex - Citation Depth:
    fetch_citation_all_depth: Iterates over an updating list of papers to process, yielding the
        response from collect papers for each paper in the list. As papers are collected, they are
        added to the set, while new, one-level deeper papers, are added to the list.
    fetch_citation_specific_depth: Iterates over an updating list of papers to process, yielding
        the response from collect papers for each paper in the list. As papers are collected, they
        are added to the set, while new, one-level deeper papers, are added to the list. Now
        includes level tracking to limit depth of citation tree exploration.
    _process_flatten_dict: Process and flatten, removing entries before 2021 and cleaning IDs.

OpenAlex - Biopython MeSH tagging:
    data_generator: Load a batch of data from the input dataset.
    json_loader: Load JSON data, transform it into a DataFrame, and perform data transformations.
    fetch_additional_mesh: Fetches additional mesh terms for DOIs in a DataFrame.
"""

import logging
from typing import Dict, List, Tuple, Generator
import pandas as pd
from joblib import Parallel, delayed
from kedro.io import AbstractDataset
from Bio import Entrez
from ..data_collection_oa.nodes import collect_papers  # pylint: disable=E0402

logger = logging.getLogger(__name__)


def fetch_citation_all_depth(
    seed_paper: str, api_config: Dict[str, str], filter_config: str
) -> Generator[Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]], None, None]:
    """
    Iterates over an updating list of papers to process, yielding the response from collect
    papers for each paper in the list. As papers are collected, they are added to the set, while
    new, one-level deeper papers, are added to the list.

    Args:
        seed_paper (str): The seed work ID paper to start from, ie. AlphaFold's.
        api_config (Dict[str, str]): The API configuration.
        filter_config (str): The filter to apply when fetching papers.

    Yields:
        Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]]: A tuple containing the edges and
        the papers.
    """
    processed_paper_ids = set()
    papers_to_process = {seed_paper}  # Use a set for uniqueness and efficient look-up

    while papers_to_process:
        logger.info("Processing %s papers", len(papers_to_process))
        current_batch = set()

        # take up to 200 papers for processing
        while papers_to_process and len(current_batch) < 200:
            current_batch.add(papers_to_process.pop())

        processed_paper_ids.update(current_batch)
        logger.info("Processing the following papers: %s", current_batch)

        # parallel fetching of papers
        child_papers = Parallel(n_jobs=6, backend="loky", verbose=10)(
            delayed(collect_papers)(
                oa_ids=paper,
                mailto=api_config["mailto"],
                perpage=api_config["perpage"],
                filter_criteria=filter_config,
                eager_loading=True,
            )
            for paper in current_batch
        )

        # flatten the list of dicts into a single dict
        child_papers_flat = {k: v for d in child_papers for k, v in d.items()}

        lengths = {key: len(value) for key, value in child_papers_flat.items()}
        logger.info("Lengths of value lists: %s", lengths)

        new_papers = set()
        for parent, children in child_papers_flat.items():
            edge_list = []

            # removing papers published before 2021-01-01
            children = [
                child
                for child in children
                if child.get("publication_date") >= "2021-01-01"
            ]

            # adding edges
            edge_list = [
                (parent, child.get("id", "").replace("https://openalex.org/", ""))
                for child in children
            ]

            # updating the list of new papers
            new_papers.update(
                [
                    clean_id
                    for child in children
                    if (
                        clean_id := child.get("id", "").replace(
                            "https://openalex.org/", ""
                        )
                    )
                    not in processed_paper_ids
                ]
            )

            edge_list_df = pd.DataFrame(edge_list, columns=["target", "source"])
            yield {parent: edge_list_df}, {parent: children}

        # extend the papers_to_process without duplicates
        papers_to_process.update(new_papers - processed_paper_ids)


def fetch_citation_to_specific_depth(
    seed_paper: str, api_config: Dict[str, str], filter_config: str, max_depth: int
):
    """
    Fetches citations to a specific depth from a seed paper.

    Args:
        seed_paper (str): The seed paper to start fetching citations from.
        api_config (Dict[str, str]): API configuration parameters.
        filter_config (str): Filter configuration for fetching papers.
        max_depth (int): The maximum depth to fetch citations to.

    Yields:
        dict: A dictionary containing the fetched papers at each level of depth.
    """
    processed_papers, papers_to_process, level = set(), {seed_paper}, 0

    while level < max_depth:
        next_level_papers = set()
        # outputs = {}

        while papers_to_process:
            current_batch = {
                papers_to_process.pop() for _ in range(min(800, len(papers_to_process)))
            }
            logger.info("Processing %d papers", len(current_batch))

            child_papers = Parallel(n_jobs=8, backend="loky", verbose=10)(
                delayed(collect_papers)(
                    oa_ids=paper,
                    mailto=api_config["mailto"],
                    perpage=api_config["perpage"],
                    filter_criteria=filter_config,
                    eager_loading=True,
                )
                for paper in current_batch
            )

            processed_papers.update(current_batch)
            child_papers_flat = {k: v for d in child_papers for k, v in d.items()}
            new_papers = _process_flatten_dict(child_papers_flat)

            next_level_papers.update(
                source["id"]
                for sources in new_papers.values()
                for source in sources
                if source["id"] not in processed_papers
            )
            # outputs.update(new_papers)
            yield {f"{level}/{key}": value for key, value in new_papers.items()}

        papers_to_process = next_level_papers
        level += 1


def _process_flatten_dict(child_papers_dict: Dict[str, any]) -> Dict[str, any]:
    """Process and flatten, removing entries before 2021 and cleaning IDs."""
    return {
        target: [
            {
                **{k: v for k, v in source.items() if k != "referenced_works"},
                "id": source.get("id", "").replace("https://openalex.org/", ""),
            }
            for source in sources
            if source.get("publication_date", "0") >= "2021-01-01"
        ]
        for target, sources in child_papers_dict.items()
        if any(
            source.get("publication_date", "0") >= "2021-01-01" for source in sources
        )
    }


# enriching the data with additional mesh tags
def data_generator(
    data: Dict[str, AbstractDataset],
    level: int,
) -> Generator[pd.DataFrame, None, None]:
    """
    Load a batch of data from the input dataset.
    """
    # select keys that start by the level
    keys = [key for key in data.keys() if key.startswith(str(level))]

    for key in keys:
        yield json_loader(data, key, level)


def json_loader(data: Dict[str, AbstractDataset], key: str, level: int) -> pd.DataFrame:
    """
    Load JSON data, transform it into a DataFrame, and perform data transformations.

    Args:
        data (Dict[str, AbstractDataset]): The input dataset.
        key (str): The key to load from the input dataset.
        level (int): The level of the data.

    Returns:
        pandas.DataFrame: The transformed DataFrame.

    """
    parent_id = key.lstrip(f"{level}/")
    # load json and keep id, doi, publication_date, mesh_terms if exist
    raw_json_data = data[key]()
    json_data = [
        {
            k: v
            for k, v in item.items()
            if k in ["id", "doi", "publication_date", "mesh_terms"]
        }
        for item in raw_json_data
    ]

    # transform to dataframe and add parent_id
    df = pd.DataFrame(json_data)
    df["parent_id"] = parent_id
    df["level"] = level

    # mesh terms is a list of dictionaries, I want to keep only a list of tuples with "descriptor_ui" and "descriptor_name" for each
    df["mesh_terms"] = df["mesh_terms"].apply(
        lambda x: [(c["descriptor_ui"], c["descriptor_name"]) for c in x] if x else None
    )

    # change doi to remove the url
    df["doi"] = df["doi"].str.replace("https://doi.org/", "")

    return df


def fetch_additional_mesh(df: pd.DataFrame) -> pd.DataFrame:
    """
    Fetches additional mesh terms for DOIs in a DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame containing DOIs and mesh terms.

    Returns:
        pd.DataFrame: The updated DataFrame with additional mesh terms.

    """
    # Filter DataFrame for rows where mesh_terms is empty
    df_no_mesh = df[df["mesh_terms"].isna()]
    logger.info("Fetching additional mesh terms for %s DOIs", len(df_no_mesh))

    # Create list of DOIs
    dois = df_no_mesh["doi"] + "[DOI]"
    # Fetch PubMed entries for each DOI
    for i, doi in enumerate(dois):
        logger.info(
            "Fetching mesh terms for DOI %s. %s/%s",
            doi,
            i + 1,
            len(dois),
        )
        stream = Entrez.esearch(db="pubmed", term=doi, retmax="1")
        record = Entrez.read(stream)

        # if I get a record, I want to extract the pubmed ID, which is the first ID in the list
        pmid = record["IdList"][0] if record["IdList"] else None

        # if pubmed_id is not None, I want to fetch the mesh terms
        if pmid:
            stream = Entrez.efetch(db="pubmed", id=pmid, retmax="250")
            record = Entrez.read(stream)

            # get the meshheadings & descriptors
            mesh_terms = record["PubmedArticle"][0]["MedlineCitation"].get(
                "MeshHeadingList", None
            )
            if mesh_terms:
                descriptors = [
                    (
                        heading["DescriptorName"].attributes["UI"],
                        str(heading["DescriptorName"]),
                    )
                    for heading in mesh_terms
                ]
            else:
                descriptors = None

            # update the dataframe with the new mesh terms for the corresponding DOI, removing the [DOI] part
            df.at[df[df["doi"] == doi.replace("[DOI]", "")].index[0], "mesh_terms"] = (
                descriptors
            )

    logger.info(
        "Still missing mesh terms for %s DOIs", len(df[df["mesh_terms"].isna()])
    )
    return df


def process_data_by_level(data: Dict[str, AbstractDataset], level: int) -> pd.DataFrame:
    """
    Process data by level and return a DataFrame with selected columns.

    Args:
        data (Dict[str, AbstractDataset]): A dictionary containing the input data.
        level (int): The level to process the data for.

    Returns:
        pd.DataFrame: A DataFrame containing the processed data with selected columns.
    """

    logger.info("Processing data for level %s", level)

    # Generate data for current level
    data_gen = data_generator(data, level)

    # Initialize an empty DataFrame for current level
    df_level = pd.DataFrame()

    # Iterate over data batches in current level
    for df_batch in data_gen:
        logger.info("Processing parent IDs: %s", df_batch["parent_id"].iloc[0])
        # Fetch additional mesh terms for current batch
        df_batch = fetch_additional_mesh(df_batch)

        # Append current batch to level DataFrame
        df_level = pd.concat([df_level, df_batch])

    return df_level[
        ["parent_id", "id", "level", "doi", "publication_date", "mesh_terms"]
    ]
