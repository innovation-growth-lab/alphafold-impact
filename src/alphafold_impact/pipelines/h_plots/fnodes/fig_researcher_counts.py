"""
This script creates the charts for the figure "Researcher Counts".
"""

import pandas as pd
import altair as alt
from ..utils import normalise_january_counts  # pylint: disable=E0402

SOURCE_LABELS = {
    "af": "AlphaFold Papers",
    "ct_ai": "AI-intensive Frontiers",
    "ct_pp": "Protein Prediction Frontiers",
    "ct_sb": "Other Struct. Biol. Frontiers",
    "other": "Other Struct. Biol. Research",
}


def create_monthly_author_counts(filtered_publications: pd.DataFrame) -> pd.DataFrame:
    """create cumulative counts of authors by month, colour by source, and
    have Adjacent and Downstream side-by-side. Counts are weighted by the
    proportion of strong users per source-and-type combination."""

    filtered_publications = filtered_publications.explode("author")
    unique_authors = filtered_publications.sort_values(
        ["strong", "weak", "publication_date"], ascending=[False, True, True]
    ).drop_duplicates(subset=["author"], keep="first")

    # Calculate strong user proportions per source-and-type combination
    # First, we need to identify the source for each author
    # Assuming there's a 'source' column in the data
    strong_proportions = (
        unique_authors.groupby(["adjacent", "downstream"])
        .agg({"strong": "sum", "weak": "sum"})  # Total strong users  # Total weak users
        .reset_index()
    )

    # Calculate proportion of strong users: strong / (strong + weak)
    strong_proportions["strong_proportion"] = strong_proportions["strong"] / (
        strong_proportions["strong"] + strong_proportions["weak"]
    )

    # Handle division by zero (set to 0 if no users)
    strong_proportions["strong_proportion"] = strong_proportions[
        "strong_proportion"
    ].fillna(0)

    # create monthly counts
    monthly_counts = (
        unique_authors.groupby(["year", "month", "adjacent", "downstream", "source"])
        .size()
        .reset_index(name="count")
    )

    # Merge with strong proportions
    monthly_counts = monthly_counts.merge(
        strong_proportions[["adjacent", "downstream", "strong_proportion"]],
        on=["adjacent", "downstream"],
        how="left",
    )

    # Apply strong user proportion weighting to counts
    monthly_counts["count"] = (
        monthly_counts["count"] * monthly_counts["strong_proportion"]
    )

    normalised_counts = normalise_january_counts(monthly_counts)

    normalised_counts["publication_date"] = pd.to_datetime(
        normalised_counts[["year", "month"]].assign(day=1)
    )

    normalised_counts["type"] = normalised_counts.apply(
        lambda x: "Adjacent" if x["adjacent"] else "Downstream", axis=1
    )

    # sort by publication date
    normalised_counts = normalised_counts.sort_values("publication_date")

    # calculate cumsums
    normalised_counts["cumsum"] = normalised_counts.groupby(["type", "source"])[
        "count"
    ].cumsum()

    # drop strong-proportion column
    normalised_counts = normalised_counts.drop(columns=["strong_proportion"])

    return normalised_counts


def create_cumul_author_charts(processed_data: pd.DataFrame) -> alt.VConcatChart:
    """Create the Altair chart for Figure Researcher Counts."""  # pylint: disable=C0301
    # Add a column for display labels
    processed_data = processed_data.copy()
    processed_data["source_label"] = processed_data["source"].map(SOURCE_LABELS)

    # Filter for Adjacent and Downstream types
    adjacent_data = processed_data[processed_data["type"] == "Adjacent"]
    downstream_data = processed_data[processed_data["type"] == "Downstream"]

    base_kwargs = dict(
        width=350,
        height=250,
    )

    x_encoding = alt.X(
        "publication_date:T",
        title="Date",
        axis=alt.Axis(format="%b %Y", labelAngle=0, labelFontSize=18),
    )
    y_encoding = alt.Y(
        "cumsum:Q",
        title="Unique Authors",
        axis=alt.Axis(
            labelFontSize=18,
            tickCount=5,
            grid=False,
        ),
    )
    color_encoding = alt.Color(
        "source_label:N",
        scale=alt.Scale(
            domain=list(SOURCE_LABELS.values()),
        ),
        legend=alt.Legend(
            title="Source",
            titleFontSize=18,
            labelFontSize=18,
            symbolSize=150,
            labelLimit=1500,
        ),
    )

    chart_adjacent = (
        alt.Chart(adjacent_data)
        .mark_line(
            point={"size": 60},
            strokeWidth=3,
            size=40,
        )
        .encode(
            x=x_encoding,
            y=y_encoding,
            color=color_encoding,
        )
        .properties(
            **base_kwargs,
            title=alt.TitleParams(
                text="(a) Cumulative count of adjacent researchers",
                anchor="start",
                fontSize=20,
            ),
        )
    )

    chart_downstream = (
        alt.Chart(downstream_data)
        .mark_line(
            point={"size": 120},
            strokeWidth=3,
            size=60,
        )
        .encode(
            x=x_encoding,
            y=y_encoding,
            color=color_encoding,
        )
        .properties(
            **base_kwargs,
            title=alt.TitleParams(
                text="(b) Cumulative count of downstream researchers",
                anchor="start",
                fontSize=20,
            ),
        )
    )

    # vconcat the two charts
    return alt.vconcat(chart_adjacent, chart_downstream).resolve_scale(y="independent")
