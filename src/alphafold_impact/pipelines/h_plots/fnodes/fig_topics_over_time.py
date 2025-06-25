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

MAIN_TOPICS = [
    "Biochemistry, Genetics and Molecular Biology",
    "Medicine",
    "Agricultural and Biological Sciences",
    "Computer Science",
    "Chemistry",
    "Immunology and Microbiology",
    "Environmental Science",
    "Materials Science",
    "Neuroscience",
    "Engineering",
    "Other",
]


def preprocess_data_topics_over_time(publications: pd.DataFrame) -> pd.DataFrame:
    """Preprocess data for Figure 2a."""
    publications = publications.copy()
    publications = publications[publications["source"] == "af"]
    publications["publication_date"] = pd.to_datetime(publications["publication_date"])
    publications["year"] = publications["publication_date"].dt.year
    publications["month"] = publications["publication_date"].dt.month

    # create boolean columns for adjacent and downstream
    publications["adjacent"] = publications["level"].eq("0")
    publications["downstream"] = publications["level"].ne("0")

    filtered_publications = publications[
        (publications["year"] >= 2021)
        & ~((publications["year"] == 2025) & (publications["month"] >= 4))
        & ~((publications["year"] == 2021) & (publications["month"] < 6))
    ]

    filtered_publications["location"] = np.select(
        [filtered_publications["adjacent"], filtered_publications["downstream"]],
        ["Adjacent", "Downstream"],
        default="unknown",
    )

    # Create a publication_date column for time series
    filtered_publications["publication_date"] = pd.to_datetime(
        filtered_publications[["year", "month"]].assign(day=1)
    )

    # Map primary fields to main topics or "Other"
    main_topics_list = MAIN_TOPICS[:-1]  # Exclude "Other" from the list
    filtered_publications["topic_group"] = filtered_publications["primary_field"].apply(
        lambda x: x if x in main_topics_list else "Other"
    )

    # group by location, topic group, and time, count the number of publications
    field_counts = (
        filtered_publications.groupby(["location", "topic_group", "publication_date"])
        .size()
        .reset_index(name="count")
    )

    # Rename topic_group back to primary_field for chart compatibility
    field_counts = field_counts.rename(columns={"topic_group": "primary_field"})

    # transform to within-location-time shares
    field_counts["share"] = field_counts.groupby(["location", "publication_date"])[
        "count"
    ].transform(lambda x: x / x.sum())

    return field_counts


def create_chart_topics_over_time(data: pd.DataFrame, title: str) -> alt.Chart:
    """Create a stacked bar chart of topics over time, faceted by location."""
    chart = (
        alt.Chart(data)
        .mark_bar(
            size=25,  # Increase bar width
        )
        .properties(
            width=300,
            height=200,
        )
        .encode(
            x=alt.X(
                "publication_date:T",
                title="Date",
                axis=alt.Axis(format="%b %Y", labelAngle=0, labelFontSize=18),
            ),
            y=alt.Y(
                "share:Q",
                title="Share",
                stack="normalize",  # Stacked and normalized to 1 (proportion)
                axis=alt.Axis(
                    labelFontSize=18,
                    tickCount=5,
                    grid=False,
                ),
            ),
            color=alt.Color(
                "primary_field:N",
                title=None,
                sort="ascending",
                legend=alt.Legend(
                    title=None,
                    labelFontSize=18,
                    titleFontSize=18,
                    labelLimit=1500,
                ),
            ),
            order=alt.Order("primary_field:N", sort="ascending"),
        )
        .facet(
            column=alt.Column(
                "location:N",
                title=None,
                header=alt.Header(labelFontSize=18, titleFontSize=18),
            )
        )
        .resolve_scale(color="shared")
        .properties(
            title=alt.TitleParams(text=title, anchor="middle"),
        )
    )
    return chart
