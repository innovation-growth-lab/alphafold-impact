"""
This is a boilerplate pipeline 'data_collection_ecr'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
from typing import Dict, List
from joblib import Parallel, delayed
from itertools import chain
from ..data_collection_oa.nodes import collect_papers  # pylint: disable=E0402
from ..data_collection_sb_labs.nodes import (  # pylint: disable=E0402
    separate_ct_from_seed,
)

logger = logging.getLogger(__name__)


def _get_sb_candidate_authors(
    data: pd.DataFrame, add_institution: bool = False
) -> pd.DataFrame:
    """
    Get the candidate authors from the given DataFrame.

    Args:
        data (pd.DataFrame): The input DataFrame containing the data.

    Returns:
        pd.DataFrame: A DataFrame containing the candidate authors and their counts.
    """

    # get the authors
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
        .drop_duplicates(subset=["id", "author"])
        .reset_index(drop=True)
    )

    # drop A9999999999
    author_data = author_data[~(author_data["author"] == "A9999999999")]

    # groupby counts of author rows
    if add_institution:
        author_data = (
            author_data.groupby(["author", "institution"])
            .size()
            .reset_index(name="counts")
        )
    else:
        author_data = author_data.groupby(["author"]).size().reset_index(name="counts")

    return author_data


def get_unique_authors(
    alphafold_data: pd.DataFrame,
    baseline_data: pd.DataFrame,
    seed_data: pd.DataFrame,
    ct_data: pd.DataFrame,
):
    """
    This node takes in the dataframes for the authors of the papers
    and returns the unique authors in the dataframes.
    """
    logger.info("Getting unique authors")
    alphafold_authors = _get_sb_candidate_authors(alphafold_data, add_institution=True)
    ct_data, other_data = separate_ct_from_seed(baseline_data, seed_data, ct_data)
    ct_authors = _get_sb_candidate_authors(ct_data, add_institution=True)
    other_authors = _get_sb_candidate_authors(other_data, add_institution=True)

    alphafold_authors["source"] = "alphafold"
    ct_authors["source"] = "ct"
    other_authors["source"] = "other"

    # merge dataframes on author ids
    authors = pd.concat(
        [alphafold_authors, ct_authors, other_authors], ignore_index=True
    )

    # # if an author has at least one row with source "alphafold", assign "alphafold" to all rows for that author
    # alphafold_authors_mask = authors['author'].isin(alphafold_authors['author'])
    # authors.loc[alphafold_authors_mask, 'source'] = 'alphafold'

    # # do the same for ct (wrt other)
    # ct_authors_mask = authors['author'].isin(ct_authors['author'])
    # authors.loc[ct_authors_mask & ~alphafold_authors_mask, 'source'] = 'ct'

    # Find the institution with the highest count for each author
    top_institutions = authors.loc[authors.groupby("author")["counts"].idxmax()][
        ["author", "institution"]
    ]
    authors = authors.drop("institution", axis=1).merge(top_institutions, on="author")

    # groupby author, institution and sum
    authors = authors.groupby(["author", "institution", "source"]).sum().reset_index()

    # pivot the data
    authors = authors.pivot_table(
        index=["author", "institution"], columns="source", values="counts", fill_value=0
    ).reset_index()

    return authors

def fetch_author_outputs(
    author_ids: List[str],
    from_publication_date: str,
    api_config: Dict[str, str],
):
    author_ids = list(set(author_ids.author.tolist()))  # Remove duplicates

    # create batches of 50 authors
    author_batches = [
        "|".join(author_ids[i : i + 25]) for i in range(0, len(author_ids), 25)
    ]

    filter_batches = [[from_publication_date, batch] for batch in author_batches]

    # slice to create parallel jobs that produce slices
    slices = [filter_batches[i : i + 100] for i in range(0, len(filter_batches), 100)]

    logger.info("Fetching papers for %d author batches", len(slices))

    for i, slice_ in enumerate(slices):
        # if i < 21:
        #     logger.info("Skipping batch number %d / %d", i + 1, len(slices))
        #     continue
        logger.info("Processing batch number %d / %d", i + 1, len(slices))

        slice_results = Parallel(n_jobs=8, backend="loky", verbose=10)(
            delayed(collect_papers)(
                oa_ids=batch,
                mailto=api_config["mailto"],
                perpage=api_config["perpage"],
                filter_criteria=["from_publication_date", "author.id"],
                slice_keys=True,
                eager_loading=True,
                skip_preprocess=True,
            )
            for batch in slice_
        )

        slice_papers = [
            paper
            for sublist in [v for d in slice_results for k, v in d.items()]
            for paper in sublist
        ]

        # create a unique list of the authors in the slice
        authors = list(chain.from_iterable([x[1].split("|") for x in slice_]))

        # match papers to the authors
        slice_papers = pd.DataFrame(slice_papers)

        slice_papers["authorships"] = slice_papers["authorships"].apply(
            lambda x: (
                [
                    author["author"].get("id", "").replace("https://openalex.org/", "")
                    for author in x
                ]
                if x
                else None
            )
        )

        # create a list of topics
        slice_papers["topics"] = slice_papers["topics"].apply(
            lambda x: (
                [
                    (
                        topic.get("id", "").replace("https://openalex.org/", ""),
                        topic.get("display_name", ""),
                        topic.get("subfield", {})
                        .get("id", "")
                        .replace("https://openalex.org/", ""),
                        topic.get("subfield", {}).get("display_name", ""),
                        topic.get("field", {})
                        .get("id", "")
                        .replace("https://openalex.org/", ""),
                        topic.get("field", {}).get("display_name", ""),
                        topic.get("domain", {})
                        .get("id", "")
                        .replace("https://openalex.org/", ""),
                        topic.get("domain", {}).get("display_name", ""),
                    )
                    for topic in x
                ]
                if x
                else None
            )
        )

        # extract concepts
        slice_papers["concepts"] = slice_papers["concepts"].apply(
            lambda x: (
                [
                    (
                        concept["id"].replace("https://openalex.org/", ""),
                        concept["display_name"],
                    )
                    for concept in x
                ]
                if x
                else None
            )
        )

        # explode the authorships
        slice_papers = slice_papers.explode("authorships")

        # filter for authorships in the author list
        slice_papers = slice_papers[slice_papers["authorships"].isin(authors)]

        # normalise_citation_counts
        slice_papers = _normalise_citation_counts(slice_papers)

        # drop counts_by_year
        slice_papers.drop(columns=["counts_by_year"], inplace=True)

        yield {f"s{i}": slice_papers}


def _normalise_citation_counts(data: pd.DataFrame):

    data = data.drop_duplicates(subset=["id"])

    # explode the 'counts_by_year' column into multiple rows
    data_exploded = data.explode("counts_by_year")

    data_exploded.reset_index(drop=True, inplace=True)

    # separate the dictionary into two columns: 'year' and 'cited_by_count'
    data_exploded[["year", "cited_by_count"]] = pd.json_normalize(
        data_exploded["counts_by_year"]
    )

    # convert the 'publication_date' column to datetime and extract the year
    data_exploded["publication_date"] = pd.to_datetime(
        data_exploded["publication_date"]
    ).dt.year

    # calculate the 't' value
    data_exploded["t"] = data_exploded["year"] - data_exploded["publication_date"]

    # drop the rows with NaN values in the 't' column
    data_exploded.dropna(subset=["t"], inplace=True)

    data_exploded["t"] = data_exploded["t"].astype(int)

    # pivot the DataFrame to get the 't' columns
    t_columns = data_exploded.pivot_table(
        index="id", columns="t", values="cited_by_count", aggfunc="sum"
    )

    # remove negative columns
    t_columns = t_columns.loc[:, (t_columns.columns >= 0)]

    t_columns.columns = ["t" + str(int(i)) for i in t_columns.columns]

    # fill remaining NaN values in 't0' with 0
    t_columns["t0"] = t_columns["t0"].fillna(0)

    # merge the 't' DataFrame back to the original DataFrame
    data = data.merge(t_columns, left_on="id", right_index=True, how="left")

    return data
