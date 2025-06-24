"""
This script creates the charts for the figure "Counts and Field Shares".
"""

import logging
import pandas as pd
import numpy as np
import altair as alt
from ..utils import normalise_january_counts  # pylint: disable=E0402


TOPICS = [
    "Biochemistry, Genetics and Molecular Biology",
    "Medicine",
    "Agricultural and Biological Sciences",
    "Computer Science",
    "Immunology and Microbiology",
    "Environmental Science",
    "Materials Science",
    "Neuroscience",
    "Engineering",
    "Chemistry",
    "Biochemistry, Genetics and Molecular Biology",
    "Medicine",
    "Agricultural and Biological Sciences",
    "Computer Science",
    "Immunology and Microbiology",
    "Environmental Science",
    "Materials Science",
    "Neuroscience",
    "Engineering",
    "Chemistry",
]

proximity_domain = [
    "Strong Adjacent",
    "Weak Adjacent",
    "Strong Downstream",
    "Weak Downstream",
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


def preprocess_data_figure(publications: pd.DataFrame, source: str) -> pd.DataFrame:
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

    filtered_publications["strong"] = filtered_publications["chain_label"].str.contains(
        "strong"
    )
    filtered_publications["weak"] = filtered_publications["chain_label"].str.contains(
        "weak"
    )

    # Calculate probabilities for strong/weak from labeled data
    labeled_data = filtered_publications[
        filtered_publications["strong"] | filtered_publications["weak"]
    ]
    p_strong = labeled_data["strong"].mean()
    p_weak = labeled_data["weak"].mean()

    # Randomly assign unlabeled data based on probabilities
    unlabeled = ~(filtered_publications["strong"] | filtered_publications["weak"])
    n_unlabeled = unlabeled.sum()

    # Generate random assignments
    np.random.seed(42)
    random_assignments = np.random.choice(
        ["strong", "weak"], size=n_unlabeled, p=[p_strong, p_weak]
    )

    # Assign strength values
    filtered_publications.loc[unlabeled, "strong"] = random_assignments == "strong"
    filtered_publications.loc[unlabeled, "weak"] = random_assignments == "weak"

    filtered_publications["strength"] = filtered_publications.apply(
        lambda x: "strong" if x["strong"] else "weak", axis=1
    )

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

    return filtered_publications


def create_monthly_paper_counts(filtered_publications: pd.DataFrame) -> pd.DataFrame:
    """Create monthly paper counts."""
    # create monthly counts
    monthly_counts = (
        filtered_publications.groupby(
            ["year", "month", "adjacent", "downstream", "strength"]
        )
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


def create_topic_counts(filtered_publications: pd.DataFrame) -> pd.DataFrame:
    """Create topic counts."""

    filtered_publications = filtered_publications[
        filtered_publications["primary_field"].isin(TOPICS)
    ]

    topic_counts = (
        filtered_publications.groupby(["adjacent", "downstream", "primary_field"])
        .size()
        .reset_index()
    )

    topic_counts["type"] = topic_counts.apply(
        lambda x: "Adjacent" if x["adjacent"] else "Downstream", axis=1
    )

    topic_counts.rename(columns={0: "size"}, inplace=True)

    # create relative share (adjacent and downstream)
    topic_counts["share"] = topic_counts["size"] / topic_counts.groupby("type")[
        "size"
    ].transform("sum")

    # mult 100
    topic_counts["share"] = topic_counts["share"] * 100
    return topic_counts


def prepare_data_for_chart(monthly_counts: pd.DataFrame) -> pd.DataFrame:
    """Prepare data for chart."""
    processed_data = monthly_counts.copy()

    processed_data["Proximity"] = (
        processed_data["strength"].str.capitalize() + " " + processed_data["type"]
    )

    # Sort by date to ensure correct cumulative calculation
    processed_data = processed_data.sort_values("publication_date")

    # Calculate the cumulative count for each category
    processed_data["Cumulative Publications"] = processed_data.groupby("Proximity")[
        "count"
    ].cumsum()

    return processed_data


def create_cumul_count_chart(processed_data: pd.DataFrame) -> alt.Chart:
    """Create cumulative count chart."""
    # Create the chart
    chart = (
        alt.Chart(processed_data)
        .mark_line(
            point={"size": 120},  # Display larger points on the line marks
            strokeWidth=3,
            size=60,
        )
        .encode(
            # Map the 'publication_date' column to the x-axis
            x=alt.X(
                "publication_date:T",
                title="Date",
                # Customize the axis to show month and year
                axis=alt.Axis(format="%b %Y", labelAngle=0, labelFontSize=18),
            ),
            # Map the new 'Cumulative Publications' column to the y-axis
            y=alt.Y(
                "Cumulative Publications:Q",
                title="Total Publications",
                axis=alt.Axis(
                    labelFontSize=18,
                    tickCount=5,  # Reduce number of y-axis labels
                    grid=False,
                ),
            ),
            # Map the 'Proximity' column to the color encoding
            color=alt.Color(
                "Proximity:N",
                # Apply the custom color scheme and legend order
                scale=alt.Scale(
                    domain=proximity_domain,
                    range=["#285EA7", "#99B6D8", "#E67E22", "#F5B784"],
                ),
                legend=alt.Legend(
                    title="Proximity",
                    titleFontSize=18,
                    labelFontSize=18,
                    symbolSize=150,
                ),
            ),
        )
        .properties(
            width=625,
            height=250,
        )
    )

    chart = chart.properties(
        title=alt.TitleParams(
            text="(a) Cumulative count of publications linked to AlphaFold2",
            anchor="start",  # right align the title
            fontSize=20,
        ),
    )

    return chart


def create_topic_chart(topic_counts: pd.DataFrame) -> alt.Chart:
    """Create topic chart."""

    sort_order = (
        topic_counts[topic_counts["type"] == "Adjacent"]
        .sort_values("share", ascending=False)["primary_field"]
        .tolist()
    )
    # Create the faceted bar chart without a title
    chart = (
        alt.Chart(topic_counts)
        .mark_bar()
        .properties(
            width=300,
            height=250,
        )
        .encode(
            y=alt.Y(
                "primary_field:N",
                sort=sort_order,
                title="Primary Field",
                axis=alt.Axis(labelLimit=1500, labelFontSize=20, titlePadding=210),
            ),
            x=alt.X(
                "share:Q",
                title="Share of Papers (%)",
                axis=alt.Axis(values=list(range(0, 101, 10))),
            ),
            color=alt.value("#4c78a8"),
        )
        .facet(
            column=alt.Column(
                "type:N",
                title=None,
                header=alt.Header(labelFontSize=18, titleFontSize=18),
            )
        )
        .resolve_scale(y="shared")
    )

    # Add the title to the top of the whole faceted chart using .properties
    chart = chart.properties(
        title=alt.TitleParams(
            text="(b) Share of adjacent and downstream publications connected to AlphaFold2 by field",  # pylint: disable=C0301
            anchor="end",  # right align the title
            fontSize=20,
        ),
    )

    return chart
