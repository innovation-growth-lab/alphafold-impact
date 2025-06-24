"""
This script creates the charts for the figure "Descriptive Plots".
"""

import pandas as pd
import numpy as np
import altair as alt
from typing import List


SOURCE_LABELS = {
    "af": "AlphaFold Papers",
    "ct_ai": "AI-intensive Frontiers",
    "ct_pp": "Protein Prediction Frontiers",
    "ct_sb": "Other Struct. Biol. Frontiers",
    "other": "Other Struct. Biol. Research",
}


COLUMNS = [
    "num_publications",
    "num_uniprot_structures",
    "num_pdb_submissions",
    "num_primary_submissions",
]
VAR_LABELS = {
    "num_uniprot_structures": "Uniprot Structures",
    "num_pdb_submissions": "PDB Submissions",
    "num_primary_submissions": "Primary Submissions",
}

CHART_TITLES = [
    "(d) Uniprot Structures per Publication (4Q RA)",
    "(e) PDB Submissions per Publication (4Q RA)",
    "(f) Primary Submissions per Publication (4Q RA)",
]


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
                ],
                "stroke": [
                    "#FF5836",  # Intense Red
                    "#1F5DAD",  # Intense Blue
                    "#00B2A2",  # Intense Green
                    "#FAB61B",  # Intense Yellow
                    "#7B4FA3",  # Purple
                ],
            },
        }
    }


igl_style()
alt.theme.enable("igl_style")


def _create_quarterly_charts(structures, columns, y_titles, titles):
    charts = []
    # Get all unique quarters, sorted
    unique_quarters = sorted(structures["quarter"].unique())
    # Show every two quarters in the x-axis ticks
    x_tick_values = unique_quarters[::2]

    for column, y_label, title in zip(columns, y_titles, titles):
        scatter = (
            alt.Chart(structures)
            .mark_circle(size=60)
            .encode(
                x=alt.X(
                    "quarter:N",
                    title="Quarter",
                    axis=alt.Axis(grid=True, values=x_tick_values),
                ),
                y=alt.Y(
                    f"{column}:Q",
                    title=y_label,
                    axis=alt.Axis(
                        grid=True,
                    ),
                ),
                color=alt.Color(
                    "source:N",
                    title=None,
                    sort=[
                        "AlphaFold Papers",
                        "AI-intensive Frontiers",
                        "Protein Prediction Frontiers",
                        "Other Struct. Biol. Frontiers",
                        "Other Struct. Biol. Research",
                    ],
                ),
            )
            .properties(width=300, height=200)
        )

        line = (
            alt.Chart(structures)
            .mark_line()
            .encode(
                x=alt.X(
                    "quarter:N",
                    title="Quarter",
                    axis=alt.Axis(values=x_tick_values),
                ),
                y=alt.Y(
                    f"{column}:Q",
                    title=y_label,
                ),
                color=alt.Color(
                    "source:N",
                    title=None,
                    sort=[
                        "AlphaFold Papers",
                        "AI-intensive Frontiers",
                        "Protein Prediction Frontiers",
                        "Other Struct. Biol. Frontiers",
                        "Other Struct. Biol. Research",
                    ],
                    legend=alt.Legend(
                        title="Source",
                        titleFontSize=14,
                        labelFontSize=14,
                        labelLimit=1500,
                    ),
                ),
            )
            # .properties(width=300, height=200)
        )

        chart = (
            alt.layer(scatter, line)
            .resolve_scale(y="shared")
            .properties(title=alt.TitleParams(title, fontSize=14))
        )
        charts.append(chart)
    return charts


def add_researcher_label(data: pd.DataFrame) -> pd.DataFrame:
    """Add a researcher label to the data."""
    data = data.sort_values(by=["author", "quarter"])
    final_observation = data.groupby("author").last().reset_index()
    conditions = [
        final_observation["af"] > 0,
        final_observation["ct_ai"] > 0,
        final_observation["ct_pp"] > 0,
        final_observation["ct_sb"] > 0,
    ]

    choices = ["af", "ct_ai", "ct_pp", "ct_sb"]
    final_observation["source"] = np.select(conditions, choices, default="other")

    # merge the source classification back into the original DataFrame
    data = data.merge(final_observation[["author", "source"]], on="author", how="left")

    return data


def process_and_create_charts(
    structures: pd.DataFrame,
    columns: list,
    source_labels: dict,
    chart_titles: list,
    concat_title: str,
) -> alt.Chart:
    """Process data and create charts"""
    # create per-publication values
    for var in columns:
        structures[f"{var}_pp"] = structures[var] / structures["num_publications"]
        structures[f"{var}_rolling"] = structures.groupby("source")[
            f"{var}_pp"
        ].transform(lambda x: x.rolling(4, min_periods=1).mean())

    charts = _create_quarterly_charts(
        structures,
        [f"{var}_rolling" for var in columns if var != "num_publications"],
        [source_labels[var] for var in columns if var != "num_publications"],
        chart_titles,
    )

    hconcat = alt.hconcat(*charts).properties(
        title=alt.TitleParams(concat_title, fontSize=20, anchor="middle")
    )
    return hconcat


def create_quarterly_vals(
    data: pd.DataFrame, variables: List[str], oper: str
) -> pd.DataFrame:
    """Create quarterly values."""

    # fillna with 0
    for var in variables:
        data[var] = data[var].fillna(0)

    vals_over_time = (
        data.groupby(["quarter", "source"])
        .agg({var: oper for var in variables})
        .reset_index()
    )

    # only keep until 2024Q1
    vals_over_time = vals_over_time[
        (vals_over_time["quarter"] <= "2025Q1")
        & (vals_over_time["quarter"] >= "2020Q1")
    ]

    # delete 2021Q1 and 2021Q2 for AF
    vals_over_time = vals_over_time[
        ~((vals_over_time["quarter"] == "2020Q3") & (vals_over_time["source"] == "af"))
        & ~(
            (vals_over_time["quarter"] == "2020Q4") & (vals_over_time["source"] == "af")
        )
        & ~(
            (vals_over_time["quarter"] == "2021Q1") & (vals_over_time["source"] == "af")
        )
        & ~(
            (vals_over_time["quarter"] == "2021Q2") & (vals_over_time["source"] == "af")
        )
    ]

    vals_over_time["quarter"] = vals_over_time["quarter"].astype(str)

    # map the source labels
    vals_over_time["source"] = vals_over_time["source"].map(SOURCE_LABELS)

    # return only since 2021Q1
    vals_over_time = vals_over_time[vals_over_time["quarter"] >= "2021Q1"]

    return vals_over_time
