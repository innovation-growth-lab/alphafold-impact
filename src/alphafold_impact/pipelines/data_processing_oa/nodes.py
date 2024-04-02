"""
This is a boilerplate pipeline 'data_analysis_exploratory'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict, Callable, Generator
import pandas as pd
from kedro.io import AbstractDataset
from Bio import Entrez

logger = logging.getLogger(__name__)


# enriching the data with additional mesh tags
def _data_generator(
    data: Dict[str, AbstractDataset],
    level: int,
) -> Generator[pd.DataFrame, None, None]:
    """
    Load a batch of data from the input dataset.
    """
    # select keys that start by the level
    keys = [key for key in data.keys() if key.startswith(str(level))]

    for key in keys:
        yield _json_loader(data, key, level)


def _json_loader(
    data: Dict[str, AbstractDataset], key: str, level: int
) -> pd.DataFrame:
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
            if k in ["id", "doi", "publication_date", "mesh_terms", "cited_by_count", "authorships"]
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

    # break atuhorship nested dictionary jsons, create triplets of authorship
    df["authorships"] = df["authorships"].apply(
        lambda x: (
            [
                (
                    author["author"]["id"].replace("https://openalex.org/", ""),
                    inst["id"].replace("https://openalex.org/", ""),
                    author["author_position"],
                )
                if author["institutions"]
                else [
                    author["author"]["id"].replace("https://openalex.org/", ""),
                    "",
                    author["author_position"],
                ]
                for author in x
                for inst in author["institutions"] or [{}]
            ]
            if x
            else None
        )
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
        Entrez.email = "david.ampudia@nesta.org.uk"

        try:
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
                df.at[
                    df[df["doi"] == doi.replace("[DOI]", "")].index[0], "mesh_terms"
                ] = descriptors

        except Exception as e:  # pylint: disable=W0718
            logger.error("Error fetching mesh terms for DOI %s: %s", doi, e)

    logger.info(
        "Still missing mesh terms for %s DOIs", len(df[df["mesh_terms"].isna()])
    )
    return df


def process_data_by_level(
    data: Dict[str, AbstractDataset], level: int, extra_mesh: str = True
) -> pd.DataFrame:
    """
    Process data by level and return a DataFrame with selected columns.

    Args:
        data (Dict[str, AbstractDataset]): A dictionary containing the input data.
        level (int): The level to process the data for.
        extra_mesh (bool): Whether to enrich the data with additional mesh terms.

    Returns:
        pd.DataFrame: A DataFrame containing the processed data with selected columns.
    """

    logger.info("Processing data for level %s", level)

    # Generate data for current level
    data_gen = _data_generator(data, level)

    # Initialize an empty DataFrame for current level
    df_level = pd.DataFrame()

    # Iterate over data batches in current level
    for df_batch in data_gen:

        if extra_mesh:
            logger.info("Processing parent IDs: %s", df_batch["parent_id"].iloc[0])
            # Fetch additional mesh terms for current batch
            df_batch = fetch_additional_mesh(df_batch)

        # Append current batch to level DataFrame
        df_level = pd.concat([df_level, df_batch])

    return df_level[
        [
            "parent_id",
            "id",
            "level",
            "doi",
            "publication_date",
            "mesh_terms",
            "cited_by_count",
            "authorships"
        ]
    ]


def combine_depth_strength_level_0(
    depth_data: pd.DataFrame,
    strength_data: Dict[str, Callable[[], pd.DataFrame]],
) -> pd.DataFrame:
    """
    Combines depth data with strength data at level 0.

    Args:
        depth_data (pd.DataFrame): The depth data to be combined.
        strength_data (Dict[str, Callable[[], pd.DataFrame]]): The strength data to be combined.

    Returns:
        pd.DataFrame: The merged dataframe containing the combined data.

    """
    logger.info("Combining depth and strength data")
    depth_data["parent_doi"] = "10.1038/s41586-021-03819-2"
    strength_data_grouped = (
        strength_data.groupby(["parent_doi", "child_doi"])[
            ["influential", "intent", "context"]
        ]
        .apply(lambda x: x.to_dict(orient="records"))
        .reset_index()
        .rename(columns={0: "strength"})
    )

    # drop parent_doi
    strength_data_grouped = strength_data_grouped.drop(["parent_doi"], axis=1)

    merged_df = pd.merge(
        depth_data,
        strength_data_grouped,
        left_on="doi",
        right_on="child_doi",
        how="left",
    )
    merged_df = merged_df.drop(["child_doi"], axis=1)

    logger.info("Merging depth and strength data complete")
    return merged_df[
        [
            "parent_id",
            "id",
            "parent_doi",
            "doi",
            "level",
            "publication_date",
            "mesh_terms",
            "strength",
            "cited_by_count",
            "authorships"
        ]
    ]


def combine_depth_strength_other_levels(
    previous_depth_data: pd.DataFrame,
    depth_data: pd.DataFrame,
    strength_data: Dict[str, Callable[[], pd.DataFrame]],
) -> pd.DataFrame:
    """
    Combines depth data, strength data, and other level data into a single DataFrame.

    Args:
        previous_depth_data (pd.DataFrame): DataFrame containing previous depth data.
        depth_data (pd.DataFrame): DataFrame containing depth data.
        strength_data (Dict[str, Callable[[], pd.DataFrame]]): Dictionary of strength data.

    Returns:
        pd.DataFrame: Merged DataFrame containing combined data.

    """
    logger.info("Combining depth and strength data")
    previous_depth_data_unique = previous_depth_data.drop_duplicates(subset=["id"])
    depth_data["parent_doi"] = depth_data["parent_id"].map(
        previous_depth_data_unique.set_index("id")["doi"]
    )

    logger.info("Loading strength data. Total number of works: %s", len(strength_data))
    strength_data = pd.concat([v() for v in strength_data.values()])

    strength_data_grouped = (
        strength_data.groupby(["parent_doi", "child_doi"])[
            ["influential", "intent", "context"]
        ]
        .apply(lambda x: x.to_dict(orient="records"))
        .reset_index()
        .rename(columns={0: "strength"})
    )
    merged_df = pd.merge(
        depth_data,
        strength_data_grouped,
        left_on=["parent_doi", "doi"],
        right_on=["parent_doi", "child_doi"],
        how="left",
    )
    merged_df = merged_df.drop(["child_doi"], axis=1)
    logger.info("Merging depth and strength data complete")
    return merged_df[
        [
            "parent_id",
            "id",
            "parent_doi",
            "doi",
            "level",
            "publication_date",
            "mesh_terms",
            "strength",
            "cited_by_count",
            "authorships"
        ]
    ]


def process_subfield_data(data: Dict[str, AbstractDataset]) -> pd.DataFrame:
    """
    Process data by level without enriching the data with additional mesh terms
        and return a DataFrame with selected columns.

    Args:
        data (Dict[str, AbstractDataset]): A dictionary containing the input data.
        level (int): The level to process the data for.

    Returns:
        pd.DataFrame: A DataFrame containing the processed data with selected columns.
    """

    logger.info("Processing subfield data. Total number of laoders: %s", len(data))
    # Initialize an empty DataFrame for current level
    subfield_df = pd.DataFrame()

    for i, loader in enumerate(data.values()):
        logger.info("Processing subfield data. Loader: %s / %s", i + 1, len(data))
        raw_json_data = loader()
        json_data = [
            {
                k: v
                for k, v in item.items()
                if k
                in [
                    "id",
                    "doi",
                    "publication_date",
                    "mesh_terms",
                    "cited_by_count",
                    "authorships",
                ]
            }
            for item in raw_json_data
        ]

        # transform to dataframe and add parent_id
        df = pd.DataFrame(json_data)

        # check mesh_terms exists else skip
        if "mesh_terms" not in df.columns:
            logger.warning("Skipping subfield data processing. No mesh terms found.")
            continue

        # mesh tuples
        df["mesh_terms"] = df["mesh_terms"].apply(
            lambda x: (
                [(c["descriptor_ui"], c["descriptor_name"]) for c in x] if x else None
            )
        )

        # break atuhorship nested dictionary jsons, create triplets of authorship
        df["authorships"] = df["authorships"].apply(
            lambda x: (
                [
                    (
                        author["author"]["id"].replace("https://openalex.org/", ""),
                        inst["id"].replace("https://openalex.org/", ""),
                        author["author_position"],
                    )
                    if author["institutions"]
                    else [
                        author["author"]["id"].replace("https://openalex.org/", ""),
                        "",
                        author["author_position"],
                    ]
                    for author in x
                    for inst in author["institutions"] or [{}]
                ]
                if x
                else None
            )
        )

        # change doi to remove the url
        df["doi"] = df["doi"].str.replace("https://doi.org/", "")

        # change id to remove openalex url
        df["id"] = df["id"].str.replace("https://openalex.org/", "")

        # Append current batch to level DataFrame
        subfield_df = pd.concat([subfield_df, df])

    return subfield_df[
        ["id", "doi", "publication_date", "mesh_terms", "cited_by_count", "authorships"]
    ]
