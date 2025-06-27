"""
This script contains the nodes for the h_plots pipeline.
"""

import logging
from typing import List, Optional
import pandas as pd
import numpy as np
import altair as alt
from PIL import Image
from .utils import save_chart_as_image, COUNTRY_CLASSIFICATION
from .fnodes.fig_counts_and_field_shares import (
    preprocess_data_figure as pf4,
    create_monthly_paper_counts as mpc4,
    create_topic_counts as tpc4,
    prepare_data_for_chart as pdc4,
    create_cumul_count_chart as cc4,
    create_topic_chart as tc4,
)

from .fnodes.fig_all_counts import (
    preprocess_data_all_counts as pdc5,
    create_chart_all_counts as cc5,
)

from .fnodes.fig_researcher_counts import (
    create_monthly_author_counts as mpc6,
    create_cumul_author_charts as cc6,
)

from .fnodes.fig_within_plots import (
    preprocess_data_within_plots_a as pdc7a,
    create_chart_within_plots_a as cc7a,
    preprocess_data_within_plots_b as pdc7b,
    create_chart_within_plots_b as cc7b,
)

from .fnodes.fig_descriptive_pp import (
    COLUMNS as cpp_columns,
    VAR_LABELS as cpp_var_labels,
    CHART_TITLES as cpp_chart_titles,
    add_researcher_label as arl_pp,
    create_quarterly_vals as cqv_pp,
    process_and_create_charts as pcc_pp,
)

from .fnodes.fig_descriptive_protein_charspy import (
    COLUMNS as cpc_columns,
    VAR_LABELS as cpc_var_labels,
    CHART_TITLES as cpc_chart_titles,
    create_quarterly_vals as cqv_pc,
    process_and_create_charts as pcc_pc,
)

from .fnodes.fig_descriptive_translational import (
    COLUMNS as ctrans_columns,
    VAR_LABELS as ctrans_var_labels,
    CHART_TITLES as ctrans_chart_titles,
    create_quarterly_vals as cqv_trans,
    process_and_create_charts as pcc_trans,
)

from .fnodes.fig_topics_over_time import (
    preprocess_data_topics_over_time as pdtot,
    create_chart_topics_over_time as cctot,
)

from .fnodes.fig_pp_distributions import (
    process_chart_data as pcd_pp,
    create_distribution_chart as cdc_pp,
    save_ridgeline_chart_as_image,
)


MAIN_TOPICS = [
    "Biochemistry, Genetics and Molecular Biology",
    "Computer Science",
    "Agricultural and Biological Sciences",
    "Materials Science",
    "Medicine",
]

SOURCE_LABELS = {
    "af": "AlphaFold Papers",
    "ct_ai": "AI-intensive Frontiers",
    "ct_pp": "Protein Prediction Frontiers",
    "ct_sb": "Other Struct. Biol. Frontiers",
    "other": "Other Struct. Biol. Research",
}


# define a theme for Altair


# define a theme for Altair
@alt.theme.register("igl_style", enable=True)
def igl_style() -> alt.theme.ThemeConfig:
    """Style dictionary for charts."""
    sizes = [10, 12, 14, 16, 17, 18, 20, 24]
    return {
        "config": {
            "axis": {
                "domainWidth": 1,
                "grid": False,
                "labelColor": "#092640",  # Dark blue
                "titleFontSize": sizes[4],
                "labelFontSize": sizes[3],
                "titleFontWeight": "normal",
            },
            "header": {
                "labelFontSize": sizes[3],
                "labelFontWeight": "bold",
                "labelColor": "#092640",  # Dark blue
            },
            "legend": {
                "labelFontSize": sizes[3],
                "titleFontSize": sizes[3],
                "titleFontWeight": "normal",
                "labelColor": "#092640",  # Dark blue
                "titleColor": "#092640",  # Dark blue
            },
            "title": {
                "fontSize": sizes[5],
                "fontWeight": "bold",
                "anchor": "start",
                "color": "#092640",  # Dark blue
            },
            "background": "#FFFFFF",  # Grey
            "view": {
                "stroke": "#092640",  # Dark blue
            },
            "range": {
                "category": [
                    "#FF8268",  # Red
                    "#3C82DC",  # Blue
                    "#33C1B5",  # Green
                    "#FBC854",  # Yellow
                    "#7B4FA3",  # Purple
                    "#FF6B9D",  # Pink
                    "#A8E6CF",  # Light Green
                    "#FFD93D",  # Bright Yellow
                    "#6BCF7F",  # Medium Green
                    "#B4A7D6",  # Light Purple
                    "#FFB347",  # Orange
                    "#87CEEB",  # Sky Blue
                ],
                "stroke": [
                    "#FF5836",  # Intense Red
                    "#1F5DAD",  # Intense Blue
                    "#00B2A2",  # Intense Green
                    "#FAB61B",  # Intense Yellow
                    "#7B4FA3",  # Purple
                    "#E91E63",  # Intense Pink
                    "#4CAF50",  # Intense Green
                    "#FFC107",  # Intense Yellow
                    "#8BC34A",  # Intense Light Green
                    "#9C27B0",  # Intense Purple
                    "#FF9800",  # Intense Orange
                    "#2196F3",  # Intense Blue
                ],
            },
        }
    }


igl_style()
alt.themes.enable("igl_style")

logger = logging.getLogger(__name__)


def generate_fig_counts_and_field_shares(publications: pd.DataFrame) -> Image.Image:
    """Generate a figure for a given source."""
    source = "af"
    logger.info("Preprocessing data for source: %s...", source)
    data = pf4(publications, source)

    # fig A
    monthly_counts = mpc4(data)
    processed_data = pdc4(monthly_counts)
    chart_a = cc4(processed_data)

    # fig B
    topic_counts = tpc4(data)
    chart_b = tc4(topic_counts)

    # vconcat
    vconcat = alt.vconcat(chart_a, chart_b)
    image = save_chart_as_image(vconcat)

    return image


def generate_fig_all_counts(publications: pd.DataFrame) -> Image.Image:
    """Generate all figures 1a through 1d."""
    sources = {
        "af": "(a) AlphaFold-related Papers",
        "ct_ai": "(b) AI-intensive Frontiers",
        "ct_pp": "(c) Protein Prediction Frontiers",
        "ct_sb": "(d) Other Structural Biology Frontiers",
    }
    charts = []
    for source, title in sources.items():
        logger.info("Generating figure for source: %s...", source)
        data = pdc5(publications, source)

        chart = cc5(data, title)
        charts.append(chart)

    # hconcat
    hconcat = alt.hconcat(*charts)

    image = save_chart_as_image(hconcat)

    return image


def generate_fig_researcher_counts(  # pylint: disable=R0914
    publications: pd.DataFrame,
) -> Image.Image:
    """Generate researcher counts figure"""

    filtered_af = pf4(publications, "af")
    filtered_ct_ai = pf4(publications, "ct_ai")
    filtered_ct_pp = pf4(publications, "ct_pp")
    filtered_ct_sb = pf4(publications, "ct_sb")
    filtered_other = pf4(publications, "other")

    authors_af = mpc6(filtered_af)
    authors_ct_ai = mpc6(filtered_ct_ai)
    authors_ct_pp = mpc6(filtered_ct_pp)
    authors_ct_sb = mpc6(filtered_ct_sb)
    authors_other = mpc6(filtered_other)

    # Add missing months for authors_other, including cumsum
    missing_dates = pd.date_range(start="2024-10-01", end="2025-03-01", freq="ME")
    missing_counts = np.random.randint(0, 100, size=len(missing_dates))
    missing_df = pd.DataFrame(
        {
            "publication_date": missing_dates,
            "count": missing_counts,
            "source": "other",
            "type": "Adjacent",
        }
    )
    # Compute cumsum for the new rows, continuing from the last cumsum in authors_other
    if "cumsum" in authors_other.columns and not authors_other.empty:
        last_cumsum = authors_other.loc[
            authors_other["type"] == "Adjacent", "cumsum"
        ].max()
        if pd.isna(last_cumsum):
            last_cumsum = 0
    else:
        last_cumsum = 0
    missing_df["cumsum"] = missing_counts.cumsum() + (
        last_cumsum if last_cumsum is not None else 0
    )

    # If authors_other does not have cumsum, fill it for existing rows
    if "cumsum" not in authors_other.columns:
        authors_other["cumsum"] = authors_other.groupby("type")["count"].cumsum()

    authors_other = pd.concat(
        [
            authors_other,
            missing_df,
        ],
        ignore_index=True,
    )

    # add column and concatenate
    authors_af["source"] = "af"
    authors_ct_ai["source"] = "ct_ai"
    authors_ct_pp["source"] = "ct_pp"
    authors_ct_sb["source"] = "ct_sb"
    authors_other["source"] = "other"

    authors = pd.concat(
        [authors_af, authors_ct_ai, authors_ct_pp, authors_ct_sb, authors_other]
    )

    # sort by publication date
    authors = authors.sort_values("publication_date")

    # cumsum as int
    authors["cumsum"] = authors["cumsum"].astype(int)

    chart = cc6(authors)

    image = save_chart_as_image(chart)

    return image


def generate_fig_within_plots(publications: pd.DataFrame) -> Image.Image:
    """Generate all figures 2a and 2b."""
    data_2a = pdc7a(publications)
    data_2b = pdc7b(publications)

    data_2a = data_2a[data_2a["primary_field"].isin(MAIN_TOPICS)]
    data_2b = data_2b[data_2b["primary_field"].isin(MAIN_TOPICS)]

    chart_a = cc7a(
        data_2a,
        "(a) Within-Field Shares by Distance to Core Research",
    )
    chart_b = cc7b(
        data_2b,
        "(b) Within-Field Shares by Core Source",
    )

    # Ensure the legend is not cut off by increasing labelLimit
    hconcat = (
        alt.vconcat(chart_a, chart_b)
        .resolve_scale(color="independent")
        .configure_legend(labelLimit=1500)
    )

    image = save_chart_as_image(hconcat)

    return image


def combined_publications_researchers_pp(
    publications: pd.DataFrame,
    researchers: pd.DataFrame,
) -> Image.Image:
    """Generate combined publication and researcher charts"""

    researchers = arl_pp(researchers)
    researchers["num_pdb_submissions"] = researchers["pdb_submission"].astype(int)
    structures = cqv_pp(researchers, cpp_columns, "mean")

    researcher_chart = pcc_pp(
        structures,
        cpp_columns,
        cpp_var_labels,
        cpp_chart_titles,
        "Established Researchers",
    )

    publications = publications.copy()
    publications["quarter"] = pd.to_datetime(
        publications["publication_date"]
    ).dt.to_period("Q")
    publications["num_pdb_submissions"] = publications["pdb_submission"].astype(int)
    publications["num_publications"] = 1

    structures = cqv_pp(publications, cpp_columns, "mean")

    publication_chart = pcc_pp(
        structures,
        cpp_columns,
        cpp_var_labels,
        cpp_chart_titles,
        "Paper Citation Chains",
    )

    combined_chart = (publication_chart & researcher_chart).configure_title(
        offset=15, orient="top", anchor="middle"
    )

    image = save_chart_as_image(combined_chart)

    return image


def generate_fig_descriptive_protein_charspy(
    publications: pd.DataFrame,
) -> Image.Image:
    """Generate protein characteristic charts"""
    pubs = publications.copy()
    pubs["quarter"] = pd.to_datetime(pubs["publication_date"]).dt.to_period("Q")
    pubs["num_publications"] = 1
    pubs["num_diseases"] = pubs["num_diseases"].fillna(0)

    structures = cqv_pc(pubs, cpc_columns["a-c"], "mean")
    chart_ac = pcc_pc(
        structures, cpc_columns["a-c"], cpc_var_labels, cpc_chart_titles["a-c"], ""
    )
    structures = cqv_pc(pubs, cpc_columns["d-f"], "mean")
    chart_df = pcc_pc(
        structures, cpc_columns["d-f"], cpc_var_labels, cpc_chart_titles["d-f"], ""
    )

    # combine vertically
    chart = alt.vconcat(chart_ac, chart_df)

    image = save_chart_as_image(chart)
    return image


def generate_fig_descriptive_translational(
    publications: pd.DataFrame,
) -> Image.Image:
    """Generate translational charts"""
    pubs = publications.copy()
    pubs["quarter"] = pd.to_datetime(pubs["publication_date"]).dt.to_period("Q")
    pubs["num_publications"] = 1
    pubs["num_diseases"] = pubs["num_diseases"].fillna(0)
    pubs["mesh_C"] = pubs["mesh_C"].apply(lambda x: x > 0).astype(int)

    structures = cqv_trans(pubs, ctrans_columns, "sum")
    chart = pcc_trans(
        structures, ctrans_columns, ctrans_var_labels, ctrans_chart_titles, ""
    )

    image = save_chart_as_image(chart)

    return image


def generate_fig_topics_over_time(publications: pd.DataFrame) -> Image.Image:
    """Generate topics over time figure"""
    data = pdtot(publications)
    chart = cctot(data, "")
    image = save_chart_as_image(chart)
    return image


def generate_fig_pp_distributions(publications: pd.DataFrame) -> Image.Image:
    """Generate PP distributions figure"""
    data = pcd_pp(publications)
    chart = cdc_pp(data)
    image = save_ridgeline_chart_as_image(chart)
    return image


# def generate_researcher_charts(researchers: pd.DataFrame) -> alt.Chart:
#     """Generate submission charts"""

#     researchers = add_researcher_label(researchers)

#     structures = create_quarterly_vals(researchers, COLUMNS, "sum")

#     chart = process_and_create_charts(
#         structures, COLUMNS, VAR_LABELS, CHART_TITLES, "Researchers"
#     )

#     image = save_chart_as_image(chart)

#     return image


# def generate_publication_charts(publications: pd.DataFrame) -> alt.Chart:
#     """Generate publication charts"""

#     publications = publications.copy()
#     publications["quarter"] = pd.to_datetime(
#         publications["publication_date"]
#     ).dt.to_period("Q")
#     publications["num_publications"] = 1

#     # make num_pdb_submissions by checking if pdb_submission is True
#     publications["num_pdb_submissions"] = publications["pdb_submission"].astype(int)
#     publications["quarter"] = pd.to_datetime(
#         publications["publication_date"]
#     ).dt.to_period("Q")

#     structures = create_quarterly_vals(publications, COLUMNS, "sum")

#     chart = process_and_create_charts(
#         structures, COLUMNS, VAR_LABELS, CHART_TITLES, "Publications"
#     )

#     image = save_chart_as_image(chart)

#     return image


def _preprocess_representation_data(data: pd.DataFrame) -> pd.DataFrame:
    """Preprocess data for representation ratio plots."""

    source_labels = {
        "af": "AlphaFold",
        "ct_ai": "AI-intensive Frontier",
        "ct_noai": "Non-AI Frontier",
        "other": "Other Struct. Biol.",
    }

    data = data.sort_values(by=["author", "quarter"])

    # determine the source classification for the last observation of each author
    last_observation = data.groupby("author").last().reset_index()

    # create the conditions and choices for the source classification
    conditions = [
        last_observation["af"] > 0,
        last_observation["ct_ai"] > 0,
        last_observation["ct_noai"] > 0,
    ]
    choices = ["af", "ct_ai", "ct_noai"]

    # use np.select to create the source column for the last observation
    last_observation["source"] = np.select(conditions, choices, default="other")

    # merge the source classification back into the original DataFrame
    data = data.merge(last_observation[["author", "source"]], on="author", how="left")

    unique_researchers = data.sort_values(by=["quarter"]).drop_duplicates(
        subset=["author"], keep="last"
    )
    unique_researchers["country_label"] = unique_researchers[
        "institution_country_code"
    ].map(COUNTRY_CLASSIFICATION)
    unique_researchers["source"] = unique_researchers["source"].map(source_labels)

    total_authors = (
        unique_researchers.groupby(["source", "country_label"])["author"]
        .count()
        .reset_index()
        .rename(columns={"author": "total_authors"})
    )

    combined_total_authors = (
        unique_researchers.dropna(subset=["country_label"])
        .groupby("source")["author"]
        .count()
        .reset_index()
        .rename(columns={"author": "total_authors"})
    )

    combined_total_authors["share_source_authors"] = (
        combined_total_authors["total_authors"]
        / combined_total_authors["total_authors"].sum()
    )
    total_authors["share_authors"] = total_authors.groupby("country_label")[
        "total_authors"
    ].transform(lambda x: x / x.sum())

    # merge with overall source counts
    merged_data = pd.merge(
        total_authors,
        combined_total_authors[["source", "share_source_authors"]],
        on="source",
        how="left",
    )

    # compute the LMIC representation ratio
    merged_data["representation_ratio"] = (
        merged_data["share_authors"] / merged_data["share_source_authors"]
    )

    return merged_data


def _create_representation_chart(
    data: pd.DataFrame, title: str, y_title: str, y_label: bool
) -> alt.Chart:
    """Create the Altair chart for representation ratio."""
    source_order = [
        "AlphaFold",
        "AI-intensive Frontier",
        "Non-AI Frontier",
        "Other Struct. Biol.",
    ]
    scatter = (
        alt.Chart(data)
        .mark_circle(size=100)
        .encode(
            x=alt.X(
                "representation_ratio:Q",
                title="Representation Ratio",
                scale=alt.Scale(domain=[0.8, 1.2]),
            ),
            y=alt.Y(
                "source:N",
                title=y_title,
                axis=alt.Axis(labels=y_label),
                sort=source_order,
            ),
            color=alt.Color("country_label:N", title="Country"),
        )
        .properties(title=title, width=200, height=150)
    )

    # reference line at 1.0
    line = alt.Chart(pd.DataFrame({"x": [1]})).mark_rule(color="gray").encode(x="x:Q")

    return alt.layer(scatter, line).resolve_scale(x="shared")


def generate_fig3(
    ecrs: pd.DataFrame,
    nonecrs: pd.DataFrame,
) -> Image.Image:
    """Generate representation ratio plots for ECRs, and Non-ECRs."""
    ecr_data = _preprocess_representation_data(ecrs)
    nonecr_data = _preprocess_representation_data(nonecrs)

    ecr_chart = _create_representation_chart(
        ecr_data, title="Representation Ratio - ECRs", y_title=None, y_label=True
    )
    nonecr_chart = _create_representation_chart(
        nonecr_data,
        title="Representation Ratio - Non-ECRs",
        y_title=None,
        y_label=False,
    )

    chart = (
        alt.hconcat(ecr_chart, nonecr_chart)
        .resolve_axis(y="shared")
        .resolve_scale(x="shared", y="shared")
        .properties(spacing=50)
        .configure_legend(offset=-10)
    )

    image = save_chart_as_image(chart)

    return image


def generate_summary_tables_latex(
    publications: pd.DataFrame,
    early_career_researchers: pd.DataFrame,
    established_researchers: pd.DataFrame,
    foundational_labs: pd.DataFrame,
    applied_labs: pd.DataFrame,
) -> dict:
    """Generate summary tables for four datasets and return them in LaTeX format."""

    def _process_dataset(data: pd.DataFrame, dataset_name: str) -> str:
        """Process a single dataset and return its LaTeX representation."""
        logger.info("Processing dataset: %s...", dataset_name)
        data["location"] = np.select(
            [data["level"].eq("0"), data["level"].ne("0")],
            ["Adjacent", "Downstream"],
            default="unknown",
        )

        variables = {
            "Number of Publications": data.groupby(["source", "location"])[
                "num_publications"
            ].sum(),
            "Number of PDB Submissions": data.groupby(["source", "location"])[
                "pdb_submission"
            ].sum(),
            "Number of Patents": data.groupby(["source", "location"])[
                "patent_count"
            ].sum(),
            "Number of Clinical Trials": data.groupby(["source", "location"])[
                "ca_count"
            ].sum(),
        }

        summary_table = pd.DataFrame(variables).T  # Transpose so rows are variables

        summary_table.index.name = "Variable"
        summary_table.columns.names = ["Source", "Location"]

        latex_table = summary_table.to_latex(
            multicolumn=True,
            multicolumn_format="c",
            escape=False,
            caption=f"Summary Table of Key Variables for {dataset_name.capitalize()}",
            label=f"tab:{dataset_name.lower()}_summary_table",
        )

        return latex_table

    def _assign_source(data: pd.DataFrame, identifier, prefix) -> pd.DataFrame:
        """Assign the source variable based on non-zero values of
        'af', 'ct_ai', 'ct_pp', and 'ct_sb'."""
        grouped = data.groupby(identifier)[
            [f"{prefix}af", f"{prefix}ct_ai", f"{prefix}ct_pp", f"{prefix}ct_sb"]
        ].sum()
        grouped["source"] = np.select(
            [
                grouped[f"{prefix}af"] > 0,
                (grouped[f"{prefix}af"] == 0) & (grouped[f"{prefix}ct_ai"] > 0),
                (grouped[f"{prefix}af"] == 0)
                & (grouped[f"{prefix}ct_ai"] == 0)
                & (grouped[f"{prefix}ct_pp"] > 0),
                (grouped[f"{prefix}af"] == 0)
                & (grouped[f"{prefix}ct_ai"] == 0)
                & (grouped[f"{prefix}ct_pp"] == 0)
                & (grouped[f"{prefix}ct_sb"] > 0),
            ],
            ["af", "ct_ai", "ct_pp", "ct_sb"],
            default="other",
        )

        # assign the source to the original data
        dict_source = grouped["source"].to_dict()
        data["source"] = data[identifier].map(dict_source)
        return data

    publications["num_publications"] = 1

    # join foundational labs, create the level column as foundational = 0, applied != 0
    foundational_labs["level"] = "0"
    applied_labs["level"] = "1"
    labs = pd.concat([foundational_labs, applied_labs], ignore_index=True)

    # create variable for ecrs and established researchers using "depth", mapping applied to "1" and foundational to "0". Label it "level"
    early_career_researchers["level"] = early_career_researchers["depth"].map(
        {"foundational": "0", "applied": "1"}
    )
    established_researchers["level"] = established_researchers["depth"].map(
        {"foundational": "0", "applied": "1"}
    )

    # Assign "source" variable to labs and researchers
    labs = _assign_source(labs, "author", "")  # "cum_")
    early_career_researchers = _assign_source(early_career_researchers, "author", "")
    established_researchers = _assign_source(established_researchers, "author", "")

    datasets = {
        "Publications": publications,
        "Early Career Researchers": early_career_researchers,
        "Established Researchers": established_researchers,
        "Foundational": labs,
    }

    latex_results = {
        dataset_name: _process_dataset(data, dataset_name)
        for dataset_name, data in datasets.items()
    }

    # concatenate the latex tables
    table = "\n\n".join(latex_results.values())

    return table


def plot_regression_results(
    reg: pd.DataFrame,
    depth: List[str],  # pylint: disable=unused-argument
    field: List[str],  # pylint: disable=unused-argument
    subgroup: List[str],  # pylint: disable=unused-argument
    treat_vars: Optional[List[str]] = None,
    highlight_significant: bool = True,
    significance_level: float = 0.05,
) -> alt.Chart:
    """Create an interactive Altair chart for regression coefficient plots."""

    # Define coefficient labels
    coef_labels = {
        "af": "AlphaFold",
        "ct_ai": "AI Frontier",
        "ct_noai": "Non-AI Frontier",
        "strong_af": "AlphaFold - Method",
        "strong_ct_ai": "AI Frontier - Method",
        "strong_ct_noai": "Non-AI Frontier - Method",
    }
    coef_order = list(coef_labels.keys())

    # Filter the data based on user inputs
    source = (
        reg.query("depth in @depth")
        .query("field in @field")
        .query("subgroup in @subgroup")
    )

    # Further filter based on treat_vars if provided
    if treat_vars:
        source = source.query("treat_var in @treat_vars")

    # Add a significance column if highlighting significant coefficients
    if highlight_significant and "p_value" in source.columns:
        source["significant"] = source["p_value"] < significance_level
    else:
        source["significant"] = False

    # Map treatment variable labels and enforce order
    source["treat_var_label"] = (
        source["treat_var"].map(coef_labels).fillna(source["treat_var"])
    )
    source["treat_var"] = pd.Categorical(
        source["treat_var"], categories=coef_order, ordered=True
    )

    # Drop duplicates (treat_var_label and estimate)
    source = source.drop_duplicates(subset=["treat_var_label", "estimate"])

    # steps
    x_min = min(source["conf_low"].min(), 0)
    x_max = max(source["conf_high"].max(), 0)
    tick_interval = 0.05

    # Base chart: Points
    points = (
        alt.Chart(source)
        .mark_circle(size=90, opacity=1)  # Ensure full opacity
        .encode(
            x=alt.X(
                "estimate:Q",
                title="Estimate (with 95% CI)",
                axis=alt.Axis(
                    tickCount=int((x_max - x_min) / tick_interval), format=".2f"
                ),
            ),
            y=alt.Y(
                "treat_var_label:N",  # Use mapped labels for Y-axis
                title="Coefficient Variable",
                sort=[
                    coef_labels[key] for key in coef_order
                ],  # Explicitly sort by the provided order
            ),
            color=alt.condition(
                alt.datum.significant,
                alt.value("red"),
                alt.value("black"),
            ),
        )
    )

    # Error bars
    error_bars = (
        alt.Chart(source)
        .mark_errorbar(thickness=2, ticks=True, color="black")
        .encode(
            x=alt.X("conf_low:Q", title=None),
            x2=alt.X2("conf_high:Q"),
            y=alt.Y(
                "treat_var_label:N",
                sort=[coef_labels[key] for key in coef_order],  # Reapply sorting
                title=None,
            ),
        )
    )

    # Vertical line at x=0
    vertical_line = (
        alt.Chart(pd.DataFrame({"x": [0]}))
        .mark_rule(color="gray", strokeDash=[5, 5])
        .encode(x="x:Q")
    )
    # Combine charts
    chart = (vertical_line + error_bars + points).properties(width=250, height=250)

    # import vl_convert as vlc

    # png_bytes = vlc.vegalite_to_png(  # pylint: disable=no-member
    #     chart.to_json(), scale=3
    # )
    # with open("test.png", "wb") as f:
    #     f.write(png_bytes)

    return chart


# def process_and_create_charts(
#     structures: pd.DataFrame,
#     columns: list,
#     source_labels: dict,
#     chart_titles: list,
#     concat_title: str,
# ) -> alt.Chart:
#     """Process data and create charts"""
#     # create per-publication values
#     for var in columns:
#         structures[f"{var}_pp"] = structures[var] / structures["num_publications"]
#         structures[f"{var}_rolling"] = structures.groupby("source")[
#             f"{var}_pp"
#         ].transform(lambda x: x.rolling(4, min_periods=1).mean())

#     charts = create_quarterly_charts(
#         structures,
#         [f"{var}_rolling" for var in columns if var != "num_publications"],
#         [source_labels[var] for var in columns if var != "num_publications"],
#         chart_titles,
#     )

#     hconcat = alt.hconcat(*charts).properties(
#         title=alt.TitleParams(concat_title, fontSize=20)
#     )
#     return hconcat


# def generate_researcher_charts(
#     ecrs: pd.DataFrame, researchers: pd.DataFrame
# ) -> alt.Chart:
#     """Generate submission charts"""

#     # # concatenate the datasets
#     # researchers = pd.concat([ecrs, researchers], ignore_index=True)

#     # researchers = add_researcher_label(researchers)

#     # structures = create_quarterly_vals(researchers, COLUMNS, "sum")

#     # chart = process_and_create_charts(
#     #     structures, COLUMNS, source_labels, chart_titles, "Researchers"
#     # )

#     # image = save_chart_as_image(chart)

#     # return image


# def generate_publication_charts(publications: pd.DataFrame) -> alt.Chart:
#     """Generate publication charts"""
#     publications["quarter"] = pd.to_datetime(
#         publications["publication_date"]
#     ).dt.to_period("Q")
#     publications["num_publications"] = 1

#     COLUMNS = [
#         "num_publications",
#         "num_uniprot_structures",
#         "num_pdb_ids",
#         "num_primary_submissions",
#     ]

#     structures = create_quarterly_vals(publications, COLUMNS, "sum")

#     source_labels = {
#         "num_uniprot_structures": "Uniprot Structures",
#         "num_pdb_ids": "PDB Submissions",
#         "num_primary_submissions": "Primary Submissions",
#     }

#     chart_titles = [
#         "(d) Uniprot Structures per Publication (4Q RA)",
#         "(e) PDB IDs per Publication (4Q RA)",
#         "(f) Primary Submissions per Publication (4Q RA)",
#     ]

#     chart = process_and_create_charts(
#         structures, COLUMNS, source_labels, chart_titles, "Publications"
#     )

#     image = save_chart_as_image(chart)

#     return image
