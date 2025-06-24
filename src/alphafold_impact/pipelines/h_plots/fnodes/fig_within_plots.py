import pandas as pd
import numpy as np
import altair as alt


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

MAIN_TOPICS = [
    "Biochemistry, Genetics and Molecular Biology",
    "Computer Science",
    "Agricultural and Biological Sciences",
    "Materials Science",
    "Medicine",
]


def preprocess_data_within_plots_a(publications: pd.DataFrame) -> pd.DataFrame:
    """Preprocess data for Figure 2a."""
    publications["publication_date"] = pd.to_datetime(publications["publication_date"])
    publications["year"] = publications["publication_date"].dt.year
    publications["month"] = publications["publication_date"].dt.month

    # create boolean columns for adjacent and downstream
    publications["adjacent"] = publications["level"].eq("0")
    publications["downstream"] = publications["level"].ne("0")

    filtered_publications = publications[
        (publications["year"] >= 2021)
        & ~((publications["year"] == 2025) & (publications["month"] >= 4))
        & ~((publications["year"] == 2021) & (publications["month"] < 5))
    ]

    filtered_publications["location"] = np.select(
        [filtered_publications["adjacent"], filtered_publications["downstream"]],
        ["Adjacent", "Downstream"],
        default="unknown",
    )

    # determine the number of adjacent publications
    num_adjacent = filtered_publications[
        filtered_publications["location"] == "Adjacent"
    ].shape[0]

    # sample the same number of downstream publications
    downstream_sample = filtered_publications[
        filtered_publications["location"] == "Downstream"
    ].sample(n=num_adjacent, random_state=1)

    balanced_publications = pd.concat(
        [
            filtered_publications[filtered_publications["location"] == "Adjacent"],
            downstream_sample,
        ]
    )

    # group by location and primary field, count the number of publications
    field_counts = (
        balanced_publications.groupby(["location", "primary_field"])
        .size()
        .reset_index(name="count")
    )

    # transform to within-field shares
    field_counts["share"] = field_counts.groupby("primary_field")["count"].transform(
        lambda x: x / x.sum()
    )

    return field_counts


def preprocess_data_within_plots_b(publications: pd.DataFrame) -> pd.DataFrame:
    """Preprocess data for Figure 2b."""
    publications["publication_date"] = pd.to_datetime(publications["publication_date"])
    publications["year"] = publications["publication_date"].dt.year
    publications["month"] = publications["publication_date"].dt.month

    # create boolean columns for adjacent and downstream
    publications["adjacent"] = publications["level"].eq("0")
    publications["downstream"] = publications["level"].ne("0")

    filtered_publications = publications[
        (publications["year"] >= 2021)
        & ~((publications["year"] == 2025) & (publications["month"] >= 4))
        & ~((publications["year"] == 2021) & (publications["month"] < 7))
    ]

    # create a random sample of 50_000 publications of each source
    sample = (
        filtered_publications.groupby("source")
        .apply(lambda x: x.sample(50_000, random_state=42))
        .reset_index(drop=True)
    )

    # group by source and primary field, count the number of publications
    field_counts = (
        sample.groupby(["source", "primary_field"]).size().reset_index(name="count")
    )

    # transform to within-field shares
    field_counts["share"] = field_counts.groupby("primary_field")["count"].transform(
        lambda x: x / x.sum()
    )

    # rename sources for better readability
    field_counts["source"] = field_counts["source"].replace(
        {
            "af": "A. AlphaFold Papers",
            "ct_ai": "B. AI-intensive Frontiers",
            "ct_pp": "C. Protein Prediction Frontiers",
            "ct_sb": "D. Other Struct. Biol. Frontiers",
            "other": "E. Other Struct. Biol. Research",
        }
    )

    return field_counts


def create_chart_within_plots_a(data: pd.DataFrame, title: str) -> alt.Chart:
    """Create the Altair chart for Figure 2."""
    return (
        alt.Chart(data)
        .mark_bar()
        .encode(
            x=alt.X("share:Q", title="Share"),
            y=alt.Y(
                "primary_field:N",
                title="Primary Field",
                sort=MAIN_TOPICS,
                axis=alt.Axis(title=None, labelLimit=400),
            ),
            color=alt.Color(
                "location:N",
                title=None,
                sort=["Adjacent", "Downstream"],
            ),
        )
        .properties(
            title=alt.TitleParams(text=title, anchor="middle"), width=300, height=250
        )
    )


def create_chart_within_plots_b(data: pd.DataFrame, title: str) -> alt.Chart:
    """Create the Altair chart for Figure 2."""
    return (
        alt.Chart(data)
        .mark_bar()
        .encode(
            x=alt.X("share:Q", title="Share"),
            y=alt.Y(
                "primary_field:N",
                title="Primary Field",
                sort=MAIN_TOPICS,
                axis=alt.Axis(title=None, labelLimit=400),
            ),
            color=alt.Color(
                "source:N",
                title=None,
                sort=[
                    "AlphaFold",
                    "AI-intensive Frontier",
                    "Non-AI Frontier",
                    "Other Struct. Biol",
                ],
            ),
        )
        .properties(
            title=alt.TitleParams(text=title, anchor="middle"), width=300, height=250
        )
    )
