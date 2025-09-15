"""
This is a boilerplate pipeline 'data_collection_other_labs'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict, List, Tuple, Generator
from kedro.io import AbstractDataset
import pandas as pd

from ..f_data_collection_foundational_labs.utils import (  # pylint: disable=E0402
    compute_avg_citation_count,
    compute_sample_publication_count,
    compute_apf,
    explode_author_data,
)

logger = logging.getLogger(__name__)


def _get_applied_candidate_authors(data: pd.DataFrame) -> List[Tuple[str, str]]:
    """
    Get the last authors of the given papers.

    Args:
        data (pd.DataFrame): The input DataFrame containing the papers.

    Returns:
        List[Tuple[str, str]]: The list of last authors.
    """
    # filter data to only consider levels 1 and 2
    applied_ids = data[data["level"].astype(str).isin(["1", "2"])]["id"].unique()

    # Identify the 'id' values in 'data' where 'level' is in the other levels ("0", "seed")
    sb_ids = data[~data["level"].astype(str).isin(["1", "2"])]["id"].unique()

    # Select the rows in 'applied_ids' that are not in 'other_ids'
    data = data[data["id"].isin(applied_ids) & ~data["id"].isin(sb_ids)]

    # get the author whose third element is the last author (parse pipe-delimited authorships)
    author_data = (
        data.drop_duplicates(subset=["id"])
        .assign(
            authorships_parsed=lambda x: x["authorships"].apply(
                lambda y: [auth.split(',') for auth in y.split('|')] if isinstance(y, str) and y != '' else []
            )
        )
        .explode("authorships_parsed")
        .assign(
            author=lambda x: x["authorships_parsed"].apply(
                lambda y: y[0] if y and len(y) > 0 else None
            ),
            institution=lambda x: x["authorships_parsed"].apply(
                lambda y: y[1] if y and len(y) > 1 else None
            ),
            position=lambda x: x["authorships_parsed"].apply(
                lambda y: y[2] if y and len(y) > 2 else None
            ),
        )
        .dropna(subset=["author"])
        .drop_duplicates(subset=["id", "author", "institution", "position"])
        .reset_index(drop=True)
    )

    # drop A9999999999
    author_data = author_data[~(author_data["author"] == "A9999999999")]

    # drop empty string institution
    author_data = author_data[~(author_data["institution"] == "")]

    # calculate counts for each author, institution pair
    counts = (
        author_data[author_data["position"] == "last"]
        .groupby(["author", "institution"])["id"]
        .count()
    )

    # determine the quantile threshold
    threshold = counts.quantile(0.75)

    # filter pairs based on the threshold
    alphafold_authors = counts[counts >= threshold]

    return list(alphafold_authors.index)


def get_ai_noai_div(ct_data, seed_papers):
    """
    Get AI and non-AI papers division based on the given ct_data and seed_papers.

    Args:
        ct_data (pd.DataFrame): The DataFrame containing the citation data.
        seed_papers (pd.DataFrame): The DataFrame containing the seed papers.

    Returns:
        tuple: A tuple containing two DataFrames - ct_ai_papers and ct_noai_papers.
            - ct_ai_papers: The DataFrame containing AI papers.
            - ct_pp_papers: The DataFrame containing PP papers.
            - ct_sb_papers: The DataFrame containing SB papers.
    """

    def create_df(ct_data, ct_data_previous, level):
        """
        Create a DataFrame based on the given ct_data, ct_data_previous, and level.

        Args:
            ct_data (pd.DataFrame): The DataFrame containing the citation data.
            ct_data_previous (pd.DataFrame): The DataFrame containing the previous level.
            level (int): The level of the DataFrame.

        Returns:
            pd.DataFrame: The DataFrame filtered based on the level and parent_id.
        """
        # Filter ct_data based on level and parent_id
        ct_level_df = ct_data[
            (ct_data["level"] == str(level))
            & (ct_data["parent_id"].isin(ct_data_previous["id"]))
        ]

        return ct_level_df

    ct_level_0_ai_papers = ct_data[
        ct_data["parent_id"].isin(
            seed_papers[seed_papers["label"].str.contains("ai")]["parent_id"]
        )
    ]

    # create ct_level_0_noai_papers
    ct_level_0_pp_papers = ct_data[
        ct_data["parent_id"].isin(
            seed_papers[~seed_papers["label"].str.contains("ai")]["parent_id"]
        )
    ]

    # create ct_level_0_sb_papers
    ct_level_0_sb_papers = ct_data[
        ct_data["parent_id"].isin(
            seed_papers[seed_papers["label"].str.contains("sb")]["parent_id"]
        )
    ]

    ct_level_1_ai_papers = create_df(ct_data, ct_level_0_ai_papers, 1)
    ct_level_1_pp_papers = create_df(ct_data, ct_level_0_pp_papers, 1)
    ct_level_1_sb_papers = create_df(ct_data, ct_level_0_sb_papers, 1)

    ct_level_2_ai_papers = create_df(ct_data, ct_level_1_ai_papers, 2)
    ct_level_2_pp_papers = create_df(ct_data, ct_level_1_pp_papers, 2)
    ct_level_2_sb_papers = create_df(ct_data, ct_level_1_sb_papers, 2)

    # Concatenate the "ai", "pp", and "sb" dataframes
    ct_ai_papers = pd.concat(
        [ct_level_0_ai_papers, ct_level_1_ai_papers, ct_level_2_ai_papers]
    )
    ct_pp_papers = pd.concat(
        [ct_level_0_pp_papers, ct_level_1_pp_papers, ct_level_2_pp_papers]
    )
    ct_sb_papers = pd.concat(
        [ct_level_0_sb_papers, ct_level_1_sb_papers, ct_level_2_sb_papers]
    )

    # Drop duplicates
    ct_ai_papers = ct_ai_papers.drop_duplicates(subset=["parent_id", "id", "level"])
    ct_pp_papers = ct_pp_papers.drop_duplicates(subset=["parent_id", "id", "level"])
    ct_sb_papers = ct_sb_papers.drop_duplicates(subset=["parent_id", "id", "level"])

    return ct_ai_papers, ct_pp_papers, ct_sb_papers


def get_candidate_authors(
    alphafold_data: pd.DataFrame,
    ct_data: pd.DataFrame,
    other_data: pd.DataFrame,
    seed_papers: pd.DataFrame,
):
    """
    Get a list of candidate authors from different data sources.

    Args:
        alphafold_data (pd.DataFrame): Dataframe containing AlphaFold data.
        ct_data (pd.DataFrame): Dataframe containing CT data.
        other_data (pd.DataFrame): Dataframe containing other data.

    Returns:
        list: A list of candidate authors from different data sources.
    """
    # get the last authors from the AlphaFold data
    alphafold_authors = _get_applied_candidate_authors(alphafold_data)

    # get the last authors from the CT data
    ct_authors = _get_applied_candidate_authors(ct_data)

    # create ct_level_0_ai_papers
    ct_ai_papers, ct_pp_papers, ct_sb_papers = get_ai_noai_div(ct_data, seed_papers)

    ct_ai_authors = _get_applied_candidate_authors(ct_ai_papers)
    ct_pp_authors = _get_applied_candidate_authors(ct_pp_papers)
    ct_sb_authors = _get_applied_candidate_authors(ct_sb_papers)

    # get the last authors from the other SB data
    other_authors = _get_applied_candidate_authors(other_data)

    # create a unique list of authors
    authors = list(set(alphafold_authors + ct_authors + other_authors))

    return (
        authors,
        alphafold_authors,
        ct_authors,
        ct_ai_authors,
        ct_pp_authors,
        ct_sb_authors,
        other_authors,
    )


def create_candidates_map(
    alphafold_authors: List,
    ct_ai_authors: List,
    ct_pp_authors: List,
    ct_sb_authors: List,
    other_authors: List,
) -> pd.DataFrame:
    """
    Create a candidates map by combining the author data from different sources.

    Args:
        alphafold_authors (List): A list of authors from the AlphaFold source.
        ct_authors (List): A list of authors from the CT source.
        other_authors (List): A list of authors from other sources.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the combined author data.

    """

    # create pandas dataframes for alphafold_authors, ct_authors, other_authors
    alphafold_df = pd.DataFrame(alphafold_authors, columns=["author", "institution"])
    alphafold_df["seed"] = "alphafold"

    ct_ai_df = pd.DataFrame(ct_ai_authors, columns=["author", "institution"])
    ct_ai_df["seed"] = "ct_ai"

    ct_pp_df = pd.DataFrame(ct_pp_authors, columns=["author", "institution"])
    ct_pp_df["seed"] = "ct_pp"

    ct_sb_df = pd.DataFrame(ct_sb_authors, columns=["author", "institution"])
    ct_sb_df["seed"] = "ct_sb"

    other_df = pd.DataFrame(other_authors, columns=["author", "institution"])
    other_df["seed"] = "other"

    # combine the candidate dataframes
    candidate_data = pd.concat([alphafold_df, ct_ai_df, ct_pp_df, ct_sb_df, other_df])

    # check for duplicates in candidate_data
    # candidate_data = candidate_data.drop_duplicates(subset=["author", "institution"])

    candidate_data["seed"] = candidate_data.groupby(["author", "institution"])[
        "seed"
    ].transform(lambda x: ",".join(x.unique()))
    candidate_data = candidate_data.drop_duplicates(subset=["author", "institution"])

    # ad-hoc fixes
    # if seed contains alphafold, go to alphafold
    candidate_data["seed"] = candidate_data["seed"].apply(
        lambda x: "alphafold" if "alphafold" in x else x
    )

    # if seed contains ct_ai, go to ct_ai
    candidate_data["seed"] = candidate_data["seed"].apply(
        lambda x: "ct_ai" if "ct_ai" in x else x
    )

    # if seed contains ct_pp, go to ct_pp
    candidate_data["seed"] = candidate_data["seed"].apply(
        lambda x: "ct_pp" if "ct_pp" in x else x
    )

    # if seed contains ct_sb, go to ct_sb
    candidate_data["seed"] = candidate_data["seed"].apply(
        lambda x: "ct_sb" if "ct_sb" in x else x
    )

    return candidate_data


def calculate_lab_determinants(
    dict_loader: AbstractDataset,
    alphafold_authors: List,
    ct_ai_authors: List,
    ct_pp_authors: List,
    ct_sb_authors: List,
    other_authors: List,
) -> Generator[Dict[str, pd.DataFrame], None, None]:
    """
    This node calculates the lab determinants for each lab based on the given data.

    Parameters:
        dict_loader (AbstractDataset): The input data loader containing the lab data.
        alphafold_authors (List): The list of AlphaFold authors.
        ct_authors (List): The list of CT authors.
        other_authors (List): The list of other authors.

    Yields:
        Generator[Dict[str, pd.DataFrame], None, None]: A dictionary containing the lab
            determinants for each lab.

    """
    for i, (key, loader) in enumerate(dict_loader.items()):

        logger.info(
            "Processing data from %s. Number %s out of %s", key, i + 1, len(dict_loader)
        )

        article_list = loader()

        json_data = [
            {
                k: v
                for k, v in item.items()
                if k
                in [
                    "id",
                    "doi",
                    "publication_date",
                    "cited_by_count",
                    "authorships",
                ]
            }
            for item in article_list
        ]

        # transform to dataframe and add parent_id
        data = pd.DataFrame(json_data)

        # break authorship nested dictionary jsons, create pipe-delimited string of authorship triplets
        data["authorships"] = data["authorships"].apply(
            lambda x: (
                "|".join([
                    f"{author['author']['id'].replace('https://openalex.org/', '')},{inst['id'].replace('https://openalex.org/', '') if inst else ''},{author['author_position']}"
                    for author in x
                    for inst in author["institutions"] or [{}]
                ]) if x else None
            )
        )

        # for each dictionary, only preserve
        author_data = explode_author_data(data)

        # get mapping between candidates and chain seeds
        candidate_map = create_candidates_map(
            alphafold_authors,
            ct_ai_authors,
            ct_pp_authors,
            ct_sb_authors,
            other_authors,
        )

        # keep if author, institution pair is in candidate_map author institution columns
        author_data = author_data.merge(
            candidate_map[["author", "institution", "seed"]],
            on=["author", "institution"],
            how="inner",
        )

        # lets drop authors whose row count is less than 10 (reduces chances of reassessment)
        author_data = author_data.groupby(["author", "institution"]).filter(
            lambda x: len(x) >= 10
        )

        # compute average annual citation count
        avg_citation_count = compute_avg_citation_count(author_data)

        # compute sample publication count
        sample_publication_count = compute_sample_publication_count(author_data)

        # compute apf
        apf_data = compute_apf(author_data)

        # merge apf and avg_citation_count
        author_data = apf_data.merge(
            avg_citation_count, on=["author", "institution", "year"], how="left"
        )

        # merge with sample publication count
        author_data = author_data.merge(
            sample_publication_count, on=["author", "institution", "year"], how="left"
        )

        logger.info("Finished processing data from %s", key)

        yield {key: author_data}
