"""
This script creates the charts for the figure "All Counts".
"""

import logging
import pandas as pd
import numpy as np
import altair as alt
from ..utils import normalise_january_counts  # pylint: disable=E0402


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
                ],
                "stroke": [
                    "#FF5836",  # Intense Red
                    "#1F5DAD",  # Intense Blue
                    "#00B2A2",  # Intense Green
                    "#FAB61B",  # Intense Yellow
                ],
            },
        }
    }


igl_style()
alt.theme.enable("igl_style")

logger = logging.getLogger(__name__)


def preprocess_data_all_counts(publications: pd.DataFrame, source: str) -> pd.DataFrame:
    """Preprocess data for a given source."""
    filtered_publications = publications[publications["source"] == source].copy()

    # convert publication_date to datetime and extract year and month
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
            & (filtered_publications["month"] < 5)
        )
    ]
    # create monthly counts
    monthly_counts = (
        filtered_publications.groupby(["year", "month", "adjacent", "downstream"])
        .size()
        .reset_index(name="count")
    )

    normalised_counts = normalise_january_counts(monthly_counts)

    normalised_counts["publication_date"] = pd.to_datetime(
        normalised_counts[["year", "month"]].assign(day=1)
    )

    normalised_counts["type"] = normalised_counts.apply(
        lambda x: "Adjacent" if x["adjacent"] else "Downstream", axis=1
    )

    return normalised_counts


def create_chart_all_counts(data: pd.DataFrame, title: str) -> alt.Chart:
    """Create the Altair chart for Figure All Counts."""  # pylint: disable=C0301

    # any data count below 0, set it to 0
    data["count"] = data["count"].clip(lower=0)

    return (
        alt.Chart(data)
        .mark_bar(width=7)
        .encode(
            x=alt.X(
                "publication_date:T",
                title="Publication Date (month)",
                axis=alt.Axis(tickCount="year"),
            ),
            y=alt.Y(
                "count:Q",
                title="Count",
            ),
            color=alt.Color(
                "type:N",
                title="Publication Type",
                legend=alt.Legend(
                    orient="top-left",
                    offset=10,
                ),
            ),
        )
        .properties(
            title={
                "text": [title],
            },
            width=200,
            height=200,
        )
    )
