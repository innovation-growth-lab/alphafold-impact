"""
This module contains the code to create the distribution charts for the PP nodes.
"""

import pandas as pd
import numpy as np
import altair as alt

alt.data_transformers.enable("vegafusion")

source_labels = {
    "af": "AlphaFold Papers",
    "ct_ai": "AI-intensive Frontiers",
    "ct_pp": "Protein Prediction Frontiers",
    "ct_sb": "Other Struct. Biol. Frontiers",
    "other": "Other Struct. Biol. Research",
}


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
alt.theme.enable("igl_style")


def process_chart_data(publications: pd.DataFrame) -> pd.DataFrame:
    """Process data for the chart."""
    filtered_publications = publications.copy()
    filtered_publications["publication_date"] = pd.to_datetime(
        filtered_publications["publication_date"]
    )
    filtered_publications["year"] = filtered_publications["publication_date"].dt.year
    filtered_publications["month"] = filtered_publications["publication_date"].dt.month

    # create boolean columns for adjacent and downstream
    filtered_publications["adjacent"] = filtered_publications["level"].eq("0")
    filtered_publications["downstream"] = filtered_publications["level"].ne("0")

    filtered_publications = filtered_publications[
        (filtered_publications["year"] >= 2021)
        & ~(
            (filtered_publications["year"] == 2025)
            & (filtered_publications["month"] >= 4)
        )
        & ~(
            (filtered_publications["year"] == 2021)
            & (filtered_publications["month"] < 6)
        )
    ]

    filtered_publications["location"] = np.select(
        [filtered_publications["adjacent"], filtered_publications["downstream"]],
        ["Adjacent", "Downstream"],
        default="unknown",
    )

    # Select the score columns and source
    score_columns = ["max_tmscore", "max_fident", "max_score"]

    # Keep only necessary columns
    data_subset = filtered_publications[["source"] + score_columns].copy()

    # Melt the data to long format
    melted_data = pd.melt(
        data_subset,
        id_vars=["source"],
        value_vars=score_columns,
        var_name="score_type",
        value_name="score_value",
    )

    # Create readable labels for the score types
    score_labels = {
        "max_tmscore": "Max Template Model Score",
        "max_fident": "Max Fractional Identity",
        "max_score": "Max RCSB PDB Shape Score",
    }
    melted_data["score_type_label"] = melted_data["score_type"].map(score_labels)

    # Remove any missing values
    melted_data = melted_data.dropna(subset=["score_value"])

    return melted_data


def create_distribution_chart(
    data: pd.DataFrame, title: str = "Score Distributions by Source"
) -> alt.Chart:
    """Create modern ridgeline plots for max_tmscore, max_fident, and max_score by source."""
    # Map source codes to human-readable labels for plotting
    data = data.copy()
    data["source_label"] = data["source"].map(source_labels).fillna(data["source"])

    # Create a chart for each score type
    charts = []
    score_types = data["score_type_label"].unique()

    step = 40
    overlap = 1

    for i, score_type in enumerate(score_types):
        score_data = data[data["score_type_label"] == score_type]

        chart = (
            alt.Chart(score_data, height=step)
            .transform_bin(
                ["bin_max", "bin_min"], "score_value", bin=alt.Bin(maxbins=50)
            )
            .transform_aggregate(
                value="count()", groupby=["source_label", "bin_min", "bin_max"]
            )
            .transform_impute(
                impute="value", groupby=["source_label"], key="bin_min", value=0
            )
            .transform_joinaggregate(max_value="max(value)", groupby=["source_label"])
            .transform_calculate(normalized_value="datum.value / datum.max_value")
            .transform_calculate(sqrt_normalized_value="sqrt(datum.normalized_value)")
            .mark_area(
                interpolate="monotone",
                fillOpacity=0.8,
                stroke="lightgray",
                strokeWidth=0.5,
            )
            .encode(
                x=alt.X("bin_min:Q")
                .bin("binned")
                .title("Score Value")
                .axis(alt.Axis(labelFontSize=14, titleFontSize=16)),
                y=alt.Y("sqrt_normalized_value:Q")
                .axis(None)
                .scale(range=[step, -step * overlap]),
                fill=alt.Fill("source_label:N")
                .legend(None)
                .sort(list(source_labels.values())),
            )
            .facet(
                row=alt.Row("source_label:N")
                .title(None)
                .header(
                    labelAngle=0,
                    labelAlign="left",
                    labelFontSize=14,
                    labelPadding=5,
                    # Only show labels on the leftmost chart
                    labels=i == 0,
                )
                .sort(list(source_labels.values()))
            )
            .properties(
                title=alt.TitleParams(
                    text=score_type, anchor="end", fontSize=16, fontWeight="bold"
                ),
                bounds="flush",
            )
        )
        charts.append(chart)

    # Combine charts horizontally
    combined_chart = (
        alt.hconcat(*charts)
        .resolve_scale(x="independent", color="independent")
        .properties(
            title=alt.TitleParams(
                text=title, anchor="start", fontSize=20, fontWeight="bold"
            )
        )
        .configure_facet(spacing=0)
        .configure_view(stroke=None)
    )

    return combined_chart
