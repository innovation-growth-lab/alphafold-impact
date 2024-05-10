"""This module contains functions to generate charts and
summary tables for time to clinical citation."""

from typing import List, Dict, Callable, Any
import pandas as pd
import altair as alt
from PIL import Image
from alphafold_impact.utils.altair import altair_to_png


def create_tcc_summary_subfields(
    oa_datasets: List[pd.DataFrame],
    tcc_datasets: List[pd.DataFrame],
    subfield_names: List[str],
    mode: str,
) -> pd.DataFrame:
    """
    Create a time to clinical citation (TCC) summary DataFrame.
    Papers that cite AF should be in 0th position in the lists.

    Args:
        oa_datasets (List[pd.DataFrame]): A list of OpenAlex papers DataFrames.
        tcc_datasets (List[pd.DataFrame]): A list of Pubmed/OpenAlex matches with TCC data.
        subfield_names (List[str]): A list of subfield names corresponding to the datasets.
        mode (str): For the subfield data 'exclude' excludes papers that cite AF,
            'include' includes papers that cite AF, 'only' uses only the papers that cite AF.

    Returns:
        pd.DataFrame: A summary TCC DataFrame.
    """
    af_oa_df = oa_datasets[0]
    af_tcc_df = tcc_datasets[0]
    af_dois = af_oa_df.doi.to_list()
    if mode == "exclude":
        oa_datasets = [af_oa_df] + [
            oa_dataset.query(f"doi not in {af_dois}") for oa_dataset in oa_datasets[1:]
        ]
        tcc_datasets = [af_tcc_df] + [
            tcc_dataset.query(f"doi not in {af_dois}")
            for tcc_dataset in tcc_datasets[1:]
        ]
    elif mode == "only":
        oa_datasets = [af_oa_df] + [
            oa_dataset.query(f"doi in {af_dois}") for oa_dataset in oa_datasets[1:]
        ]
        tcc_datasets = [af_tcc_df] + [
            tcc_dataset.query(f"doi in {af_dois}") for tcc_dataset in tcc_datasets[1:]
        ]

    summary_data = []

    for subfield, oa_df, tcc_df in zip(subfield_names, oa_datasets, tcc_datasets):
        oa_df_filtered = oa_df.query("publication_date >= '2021-07-15'")
        oa_papers = len(oa_df_filtered)
        pubmed_papers = len(tcc_df)
        clinical_citations = tcc_df["num_cited_by_clin"].sum()
        cc_per_pubmed_paper = clinical_citations / pubmed_papers
        median_days_tcc = tcc_df.query("days_to_clinical_trial.notna()")[
            "days_to_clinical_trial"
        ].median()

        summary_data.append(
            {
                "Subfield": subfield,
                "OA Papers": oa_papers,
                "Pubmed Papers": pubmed_papers,
                "Clinical Citations": clinical_citations,
                "Clinical Citations per Pubmed Paper": cc_per_pubmed_paper,
                "Median Days TCC": median_days_tcc,
            }
        )

    return pd.DataFrame(summary_data)


def create_tcc_summary_af_levels(
    af_oa_data: pd.DataFrame,
    af_tcc_data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Create a time to clinical citation (TCC) summary DataFrame for
    the papers that cite AlphaFold by citation 'level'.

    Args:
        af_oa_data (pd.DataFrame): Papers in OpenAlex that cite AF (levels 0-3).
        af_tcc_data (pd.DataFrame): Pubmed/OpenAlex matches for papers that cite AF with TCC data.

    Returns:
        pd.DataFrame: A summary TCC DataFrame for the papers that cite AlphaFold by level.
    """
    oa_datasets = [
        af_oa_data,
        af_oa_data.query("level == 0"),
        af_oa_data.query("level == 1"),
        af_oa_data.query("level == 2"),
        af_oa_data.query("level == 3"),
    ]
    tcc_datasets = [
        af_tcc_data,
        af_tcc_data.query("level == 0"),
        af_tcc_data.query("level == 1"),
        af_tcc_data.query("level == 2"),
        af_tcc_data.query("level == 3"),
    ]
    level_names = [
        "Papers that cite AF (levels 0-3)",
        "Papers that cite AF (level 0)",
        "Papers that cite AF (level 1)",
        "Papers that cite AF (level 2)",
        "Papers that cite AF (level 3)",
    ]
    summary_data = []

    for subfield, oa_df, tcc_df in zip(level_names, oa_datasets, tcc_datasets):
        oa_df_filtered = oa_df.query("publication_date >= '2021-07-15'")
        oa_papers = len(oa_df_filtered)
        pubmed_papers = len(tcc_df)
        clinical_citations = tcc_df["num_cited_by_clin"].sum()
        cc_per_pubmed_paper = clinical_citations / pubmed_papers
        median_days_tcc = tcc_df.query("days_to_clinical_trial.notna()")[
            "days_to_clinical_trial"
        ].median()

        summary_data.append(
            {
                "Subfield": subfield,
                "OA Papers": oa_papers,
                "Pubmed Papers": pubmed_papers,
                "Clinical Citations": clinical_citations,
                "Clinical Citations per Pubmed Paper": cc_per_pubmed_paper,
                "Median Days TCC": median_days_tcc,
            }
        )

    return pd.DataFrame(summary_data)


def reorder_dict_by_keys(
    original_dict: Dict[Any, Any], key_order: List[Any]
) -> Dict[Any, Any]:
    """
    Reorders a dictionary based on a given list of keys.

    Args:
        original_dict (Dict[Any, Any]): The original dictionary to reorder.
        key_order (List[Any]): The list of keys defining the new order.

    Returns:
        Dict[Any, Any]: A new dictionary reordered according to `key_order`.
    """
    return {key: original_dict[key] for key in key_order if key in original_dict}


def save_tcc_summary_subfields(
    tcc_subfield_partitions: Dict[str, Callable[[], pd.DataFrame]],
    oa_datasets: List[pd.DataFrame],
    partition_order: list,
) -> Dict[str, Callable[[], pd.DataFrame]]:
    """
    Returns a dictionary of callable functions designed to generate summary TCC DataFrames.

    Args:
        tcc_subfield_partitions (Dict[str, Callable[[], pd.DataFrame]]): A dictionary where
            each key is a string representing a subfield, and each value is a callable that returns
            a DataFrame containing data for that subfield.
        oa_datasets (List[pd.DataFrame]): A list of DataFrames containing OpenAlex data
            for various subfields.
        partition_order (List[str]): A list of strings that define the new order of keys
            for organizing the 'tcc_subfield_partitions' dictionary. This controls the order
            in which summaries are generated and processed.

    Returns:
        Dict[str, Callable[[], pd.DataFrame]]: A dictionary where each key is a string
            describing the summary table (e.g., 'tcc_summary_subfields_include_af') and
            each value is a callable that generates the corresponding summary DataFrame.
    """
    tcc_subfield_partitions = reorder_dict_by_keys(
        tcc_subfield_partitions, partition_order
    )
    tcc_datasets = []
    subfield_names = []
    for subfield_name, tcc_dataset in tcc_subfield_partitions.items():
        tcc_dataset = tcc_dataset()
        tcc_datasets.append(tcc_dataset)
        subfield_names.append(subfield_name)

    return {
        "tcc_summary_subfields_include_af": lambda: create_tcc_summary_subfields(
            oa_datasets, tcc_datasets, subfield_names, mode="include"
        ),
        "tcc_summary_subfields_exclude_af": lambda: create_tcc_summary_subfields(
            oa_datasets, tcc_datasets, subfield_names, mode="exclude"
        ),
        "tcc_summary_subfields_only_af": lambda: create_tcc_summary_subfields(
            oa_datasets, tcc_datasets, subfield_names, mode="only"
        ),
        "tcc_summary_af_levels": lambda: create_tcc_summary_af_levels(
            oa_datasets[0], tcc_datasets[0]
        ),
    }


def plot_tcc_each_quarter_subfields(
    tcc_subfield_partitions: Dict[str, Callable[[], pd.DataFrame]],
    mode: str = "exclude",
    max_quarter: str = "2023Q4",
) -> Image.Image:
    """
    Plots the median days to clinical citation against quarterly publication date,
    with lines colored by the 'subfield' column.

    Args:
        tcc_subfield_partitions (Dict[str, Callable[[], Any]]): TCC data for papers
            that cite AF and subfields.
        mode (str): 'exclude' excludes papers that cite AF from the subfield data,
            any other input includes the papers that cite Af.
        max_quarter (str): Last quarter to plot. Defaults to '2023Q4'.

    Returns:
        Image.Image: PNG chart.
    """
    af_tcc_data = (
        tcc_subfield_partitions["Papers that cite AF levels 0-3"]()
        .assign(subfield="Papers that cite AF levels 0-3")
        .dropna(subset="doi")
    )
    af_dois = af_tcc_data.doi.to_list()

    combined_data = [af_tcc_data]
    tcc_subfield_partitions.pop("Papers that cite AF levels 0-3")

    for subfield_name, tcc_dataset in tcc_subfield_partitions.items():
        tcc_dataset = tcc_dataset()
        if mode == "exclude":
            tcc_dataset = tcc_dataset.query(f"doi not in {af_dois}")
        combined_data.append(tcc_dataset.assign(subfield=subfield_name))

    combined_data = pd.concat(combined_data).query("days_to_clinical_trial.notna()")

    combined_data["publication_date"] = pd.to_datetime(
        combined_data["publication_date"]
    )
    combined_data["quarter"] = combined_data["publication_date"].dt.to_period("Q")
    max_quarter = pd.Period(max_quarter, freq="Q")
    combined_data = combined_data[combined_data["quarter"] <= max_quarter]

    grouped_df = (
        combined_data.groupby(["quarter", "subfield"])
        .agg(median_days_to_clinical_trial=("days_to_clinical_trial", "median"))
        .reset_index()
    )

    # Convert 'quarter' back to string for plotting purposes
    grouped_df["quarter"] = grouped_df["quarter"].astype(str)

    chart = (
        alt.Chart(grouped_df)
        .mark_line()
        .encode(
            x=alt.X(
                "quarter:O", title="Publication quarter", axis=alt.Axis(labelAngle=0)
            ),
            y=alt.Y(
                "median_days_to_clinical_trial:Q",
                title="Median days to clinical citation",
                scale=alt.Scale(zero=True),
            ),
            color=alt.Color(
                "subfield:N",
                title="Subfield",
                legend=alt.Legend(orient="top-right", offset=10),
            ),
        )
        .properties(
            title="Median days to clinical citation by publication quarter",
            width=800,
            height=400,
        )
    )
    return altair_to_png(chart)


def plot_tcc_each_quarter_af_levels(
    tcc_subfield_partitions: Dict[str, Callable[[], pd.DataFrame]]
) -> alt.Chart:
    """
    Plots the median days to clinical citation against quarterly publication date, with lines
    colored by the citation 'level' column.

    Args:
        tcc_subfield_partitions (Dict[str, Callable[[], Any]]): TCC data for papers that cite AF and subfields.

    Returns:
        Image.Image: PNG chart.
    """
    af_tcc_data = (
        tcc_subfield_partitions["Papers that cite AF levels 0-3"]()
        .query("level.notna()")
        .query("days_to_clinical_trial.notna()")
    )
    af_tcc_data["publication_date"] = pd.to_datetime(af_tcc_data["publication_date"])
    af_tcc_data["quarter"] = af_tcc_data["publication_date"].dt.to_period("Q")
    max_quarter = pd.Period("2023Q4", freq="Q")
    af_tcc_data = af_tcc_data[af_tcc_data["quarter"] <= max_quarter]
    grouped_df = (
        af_tcc_data.groupby(["quarter", "level"])
        .agg(median_days_to_clinical_trial=("days_to_clinical_trial", "median"))
        .reset_index()
    )

    # Convert 'quarter' back to string for plotting purposes
    grouped_df["quarter"] = grouped_df["quarter"].astype(str)

    chart = (
        alt.Chart(grouped_df)
        .mark_line()
        .encode(
            x=alt.X(
                "quarter:O", title="Publication quarter", axis=alt.Axis(labelAngle=0)
            ),
            y=alt.Y(
                "median_days_to_clinical_trial:Q",
                title="Median days to clinical citation",
            ),
            color=alt.Color(
                "level:N",
                title="Level",
                legend=alt.Legend(orient="top-right", offset=10),
            ),
        )
        .properties(
            title="Median days to clinical citation by publication quarter",
            width=800,
            height=400,
        )
    )

    return altair_to_png(chart)
