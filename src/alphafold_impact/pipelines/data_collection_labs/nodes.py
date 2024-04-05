"""This module contains the nodes used in the data collection pipeline for labs.

The nodes are used to collect papers from the OpenAIRE API, extract candidate
authors from the papers, and fetch publications for the candidate authors. 
"""

import logging
from typing import Dict, List, Tuple, Generator
from kedro.io import AbstractDataset
import pandas as pd
from joblib import Parallel, delayed
from ..data_collection_oa.nodes import collect_papers  # pylint: disable=E0402

logger = logging.getLogger(__name__)


def get_candidate_authors(
    data: pd.DataFrame, baseline: bool = False
) -> List[Tuple[str, str]]:
    """
    Get the last authors of the given papers.

    Args:
        data (pd.DataFrame): The input DataFrame containing the papers.

    Returns:
        List[Tuple[str, str]]: The list of last authors.
    """
    # threshold = 10 if baseline else 1

    # get the author whose third element is the last author
    author_data = (
        data.drop_duplicates(subset=["id"])
        .explode("authorships")
        .assign(
            author=lambda x: x["authorships"].apply(
                lambda y: y[0] if y is not None else None
            ),
            institution=lambda x: x["authorships"].apply(
                lambda y: y[1] if y is not None else None
            ),
            position=lambda x: x["authorships"].apply(
                lambda y: y[2] if y is not None else None
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

    if baseline:
        # filter to keep only relevant authors (5> as last, or first author)
        author_institution_counts = author_data.groupby(
            ["author", "institution", "position"]
        )["id"].nunique()
        filtered_pairs = author_institution_counts[
            (author_institution_counts >= 5)
            & (
                author_institution_counts.index.get_level_values("position").isin(
                    ["first", "last"]
                )
            )
        ]
        baseline_authors = [
            (author, institution)
            for author, institution, _ in list(filtered_pairs.index)
        ]

        return list(set(baseline_authors))

    # instead return all author, institution pairs when position is last
    alphafold_authors = (
        author_data[author_data["position"] == "last"]
        .groupby(["author", "institution"])["id"]
        .nunique()
    )

    return list(alphafold_authors.index)


def merge_candidate_authors(data: AbstractDataset) -> List[Tuple[str, str]]:
    """
    Merges the candidate authors from the given dataset.

    Args:
        data (AbstractDataset[Dict[str, List[str]]]): The dataset containing candidate authors.

    Returns:
        List[Tuple[str, str]]: A list of tuples representing the merged candidate authors.

    """
    candidates = [item for sublist in data.values() for item in sublist()]

    return list(set(tuple(i) for i in candidates))


def fetch_author_publications(
    author_ids: List[str],
    from_publication_date: str,
    api_config: Dict[str, str],
) -> Generator[Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]], None, None]:
    """
    Fetches publications for a list of authors based on their IDs and publication date.

    Args:
        author_ids (List[str]): List of author IDs.
        from_publication_date (str): Starting publication date in the format 'YYYY-MM-DD'.
        api_config (Dict[str, str]): API configuration dictionary containing 'mailto' and 
            'perpage' values.

    Yields:
        Tuple[Dict[str, pd.DataFrame], Dict[str, List[dict]]]: A generator that yields a 
            tuple of dictionaries. The first dictionary contains dataframes with publication 
            information, where the keys are slice IDs. The second dictionary contains lists 
            of publication dictionaries, where the keys are slice IDs.

    """
    author_ids = [author[0] for author in author_ids]

    # create batches of 50 authors
    author_batches = [
        "|".join(author_ids[i : i + 50]) for i in range(0, len(author_ids), 50)
    ]

    filter_batches = [[from_publication_date, batch] for batch in author_batches]

    # slice to create parallel jobs that produce slices
    slices = [filter_batches[i : i + 25] for i in range(0, len(filter_batches), 25)]

    logger.info("Fetching papers for %d author batches", len(slices))

    for i, slice_ in enumerate(slices):

        logger.info("Processing batch number %d / %d", i + 1, len(slices))

        slice_results = Parallel(n_jobs=6, backend="loky", verbose=10)(
            delayed(collect_papers)(
                oa_ids=batch,
                mailto=api_config["mailto"],
                perpage=api_config["perpage"],
                filter_criteria=["from_publication_date", "author.id"],
                slice_keys=True,
                eager_loading=True,
                skip_preprocess=True,
                sample_size=100 * len(batch[1].split("|")),
            )
            for batch in slice_
        )

        slice_papers = [
            paper
            for sublist in [v for d in slice_results for k, v in d.items()]
            for paper in sublist
        ]

        # for each listed dict, keep "id", "authorships", "publication_date", "cited_by_count"
        slice_papers = [
            {
                "id": paper["id"].replace("https://openalex.org/", ""),
                "authorships": paper["authorships"],
                "publication_date": paper["publication_date"],
                "cited_by_count": paper["cited_by_count"],
            }
            for paper in slice_papers
        ]

        yield {f"s{i}": slice_papers}
