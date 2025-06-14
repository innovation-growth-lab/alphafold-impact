"""
This is a boilerplate pipeline 'analysis_descriptive_translational'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
from kedro.io import AbstractDataset
from ..f_data_processing_foundational_labs.utils import (  # pylint: disable=relative-beyond-top-level
    get_intent,
    get_cum_sums,
    get_strong_cum_sums,
    preprocess_for_staggered_design,
    get_pdb_activity,
    collect_covid_references,
    calculate_field_share,
    calculate_primary_field,
    calculate_mesh_balance,
    get_ai_use,
    get_patent_citations,
    get_cc,
    get_quarterly_aggregate_outputs,
)


logger = logging.getLogger(__name__)


# def get_applied_lab_staggered_outputs(
#     data_loaders: AbstractDataset,
#     mapping_df: AbstractDataset,
#     publications_data: pd.DataFrame,
#     pdb_submissions: pd.DataFrame,
#     patents_data: pd.DataFrame,
#     icite_data: pd.DataFrame,
#     mesh_terms: pd.DataFrame,
#     institutional_data: pd.DataFrame,
# ):
#     """
#     Retrieves outputs from data loaders and performs data processing.

#     Args:
#         data_loaders (AbstractDataset): A dictionary-like object containing data loaders.
#         mapping_df (AbstractDataset): A dataset containing mapping information.

#     Returns:
#         pd.DataFrame: The processed data.

#     """

#     # mesh terms extract tree_number first character
#     mesh_terms["term_group"] = mesh_terms["tree_number"].apply(
#         lambda x: str(x)[:1] if x is not None else None
#     )

#     # create dictionary of DUI to term_group
#     mesh_terms_dict = mesh_terms.set_index("DUI")["term_group"].to_dict()

#     # create dictionary of unique pi_id, intent
#     pi_intent = get_intent(publications_data)

#     # change institutional_data columns to institution, unless column is institution
#     institutional_data.drop_duplicates(subset="author", inplace=True)
#     institutional_data.columns = [
#         "institution_" + col if col != "institution" else col
#         for col in institutional_data.columns
#     ]

#     # subset applied publications
#     applied_publications = publications_data.copy().loc[
#         (publications_data["level"] != "0") & (publications_data["level"] != "-1")
#     ]

#     # extract author id (first item in each sublist in authorships)
#     applied_authors = applied_publications[["id", "authorships", "source"]].copy()
#     applied_authors["authorships"] = publications_data["authorships"].apply(
#         lambda x: [y[0] for y in x] if x is not None else []
#     )
#     applied_authors = applied_authors.explode("authorships")
#     applied_authors.drop_duplicates(subset=["id", "authorships"], inplace=True)

#     # create a dictionary of author counts, to capture the extensive margin (see for loop below)
#     author_counts = (
#         applied_authors.groupby(["authorships", "source"])
#         .size()
#         .unstack(fill_value=0)
#         .reset_index()
#     )
#     author_counts.columns = [
#         "ext_" + col if col != "authorships" else col for col in author_counts.columns
#     ]

#     outputs = []
#     agg_outputs = []
#     for i, loader in enumerate(data_loaders.values()):
#         logger.info("Loading data batch %d / %d", i + 1, len(data_loaders))
#         data_batch = loader()

#         # drop "display_name", "authorships", "ids" if they exist
#         for col in ["display_name", "authorships", "ids"]:
#             try:
#                 data_batch = data_batch.drop(columns=col)
#             except KeyError:
#                 pass

#         # drop rows with publication_date older than 2017-06-01
#         data_batch = data_batch[data_batch["publication_date"] >= "2015-06-01"]

#         # adding primary field for paper
#         data_batch["primary_field"] = data_batch["topics"].apply(
#             lambda x: x[0][5] if x is not None and len(x) > 0 else ""
#         )

#         # transform concepts to be a list of the ids (ie first item in sublists)
#         data_batch["concepts"] = data_batch["concepts"].apply(
#             lambda x: [y[0] for y in x] if x is not None else []
#         )

#         logger.info("Merging data with mapping information")
#         mapping_df.drop_duplicates(subset="author", inplace=True)
#         data_batch = data_batch.merge(
#             mapping_df[["author", "seed"]],
#             left_on="pi_id",
#             right_on="author",
#             how="left",
#         )

#         data_batch["publication_date"] = pd.to_datetime(data_batch["publication_date"])
#         data_batch["quarter"] = data_batch["publication_date"].dt.to_period("Q")

#         logger.info("Merging data with level0 data")
#         data_processed = preprocess_for_staggered_design(
#             data_batch, applied_publications
#         )

#         logger.info("Calculate cumulative sums")
#         data_processed = get_cum_sums(data_processed, applied_publications)
#         data_processed = get_strong_cum_sums(data_processed, applied_publications)

#         logger.info("Merging data with pdb_submissions data")
#         data_processed = get_pdb_activity(data_processed, pdb_submissions)

#         logger.info("Map intent links to labs")
#         data_processed["seed_to_merge"] = data_processed["seed"].apply(
#             lambda x: (
#                 "ct"
#                 if isinstance(x, str) and "ct" in x
#                 else "af" if isinstance(x, str) and "alphafold" in x else x
#             )
#         )
#         data_processed = data_processed.merge(
#             pi_intent[["pi_id", "source", "intent"]],
#             left_on=["pi_id", "seed_to_merge"],
#             right_on=["pi_id", "source"],
#             how="left",
#         )
#         data_processed.drop(columns=["seed_to_merge", "source"], inplace=True)

#         logger.info("Map papers with protein_prediction links")
#         data_processed["protein_concept"] = data_processed["concepts"].apply(
#             lambda x: (
#                 any(element == "C18051474" or element == "C47701112" for element in x)
#             )
#         )

#         logger.info("Check if paper includes AI concept")
#         data_processed["ai_concept"] = data_processed["concepts"].apply(get_ai_use)

#         logger.info("Map papers with experimental links")
#         data_processed["experimental"] = data_processed["concepts"].apply(
#             lambda x: (
#                 any(
#                     element == "C185592680"
#                     or element == "C55493867"
#                     or element == "C12554922"
#                     or element == "C46141821"
#                     or element == "C121332964"
#                     for element in x
#                 )
#             )
#         )

#         logger.info("Collect COVID 2020 references")
#         data_processed = collect_covid_references(data_processed)

#         logger.info("Extracting fields")
#         data_processed = calculate_field_share(data_processed)

#         logger.info("Extracting primary field")
#         data_processed = calculate_primary_field(data_processed)

#         logger.info("Extracting mesh tags")
#         data_processed = calculate_mesh_balance(data_processed, mesh_terms_dict)

#         logger.info("Merging data with patents data")
#         data_processed = get_patent_citations(data_processed, patents_data)

#         logger.info("Merging data with clinical citations data")
#         data_processed = get_cc(data_processed, icite_data)

#         logger.info("Compute the extensive margin of citations to af, ct, other")
#         data_processed = data_processed.merge(
#             author_counts, left_on="pi_id", right_on="authorships", how="left"
#         )

#         logger.info("Discard columns not needed for analysis")
#         data_processed = data_processed.drop(
#             columns=["author", "topics", "mesh_terms", "concepts", "counts_by_year"]
#         )

#         logger.info("Merging data with institutional information")
#         data_processed = data_processed.merge(
#             institutional_data,
#             left_on="pi_id",
#             right_on="institution_author",
#             how="left",
#         )

#         # drop institution_author
#         data_processed = data_processed.drop(columns=["institution_author"])

#         outputs.append(data_processed)
#         quarterly = get_quarterly_aggregate_outputs(data_processed)
#         agg_outputs.append(quarterly)

#     data = pd.concat(outputs)
#     agg_data = pd.concat(agg_outputs)

#     return data, agg_data
