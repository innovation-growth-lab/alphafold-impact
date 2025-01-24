import logging
import pandas as pd
import numpy as np
import altair as alt
from PIL import Image
from .utils import normalise_january_counts, save_chart_as_image, COUNTRY_CLASSIFICATION


MAIN_TOPICS = [
    "Biochemistry, Genetics and Molecular Biology",
    "Computer Science",
    "Agricultural and Biological Sciences",
    "Materials Science",
    "Medicine",
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
                "titleFontSize": sizes[2],
                "labelFontSize": sizes[1],
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
                "fontSize": sizes[4],
                "fontWeight": "bold",
                "anchor": "start",
                "color": "#092640",  # Dark blue
            },
            "background": "#EBEBEB",  # Grey
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
alt.themes.enable("igl_style")

logger = logging.getLogger(__name__)


def _preprocess_data_fig1(publications: pd.DataFrame, source: str) -> pd.DataFrame:
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

    # create monthly counts
    monthly_counts = (
        filtered_publications.groupby(["year", "month", "adjacent", "downstream"])
        .size()
        .reset_index(name="count")
    )

    monthly_counts = monthly_counts[
        (monthly_counts["year"] >= 2021)
        & (monthly_counts["year"] != 2025)
        & ~((monthly_counts["year"] == 2024) & (monthly_counts["month"] >= 5))
        & ~((monthly_counts["year"] == 2021) & (monthly_counts["month"] == 1))
    ]

    normalised_counts = normalise_january_counts(monthly_counts)

    normalised_counts["publication_date"] = pd.to_datetime(
        normalised_counts[["year", "month"]].assign(day=1)
    )

    normalised_counts["type"] = normalised_counts.apply(
        lambda x: "Adjacent" if x["adjacent"] else "Downstream", axis=1
    )

    return normalised_counts


def _create_chart_fig1(data: pd.DataFrame, title: str) -> alt.Chart:
    """Create the Altair chart for Figure 1."""

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
            width=275,
            height=200,
        )
    )


def _fig1_process(publications: pd.DataFrame, source: str, title: str) -> Image.Image:
    """Generate a figure for a given source."""
    logger.info("Preprocessing data for source: %s...", source)
    data = _preprocess_data_fig1(publications, source)

    logger.info("Creating chart for source: %s...", source)
    chart = _create_chart_fig1(data, title)

    logger.info("Saving chart as an image for source: %s...", source)
    image = save_chart_as_image(chart)

    return chart, image


def generate_fig1(publications: pd.DataFrame) -> Image.Image:
    """Generate all figures 1a through 1d."""
    sources = {
        "af": "(a) Count of AlphaFold-related Papers",
        "ct_ai": "(b) Count of AI-intensive Frontier Papers",
        "ct_noai": "(c) Count of non-AI Frontier Papers",
        "other": "(d) Count of Other Structural Biology Papers",
    }
    figures = {}
    charts = []
    for source, title in sources.items():
        logger.info("Generating figure for source: %s...", source)
        chart, figures[source] = _fig1_process(publications, source, title)
        charts.append(chart)

    # hconcat
    hconcat1 = alt.hconcat(charts[0], charts[1])
    hconcat2 = alt.hconcat(charts[2], charts[3])
    vconcat = alt.vconcat(hconcat1, hconcat2)

    image = save_chart_as_image(vconcat)

    return image


def _preprocess_data_fig2a(publications: pd.DataFrame) -> pd.DataFrame:

    publications["location"] = np.select(
        [publications["level"].eq("0"), publications["level"].ne("0")],
        ["Adjacent", "Downstream"],
        default="unknown",
    )

    # determine the number of adjacent publications
    num_adjacent = publications[publications["location"] == "Adjacent"].shape[0]

    # sample the same number of downstream publications
    downstream_sample = publications[publications["location"] == "Downstream"].sample(
        n=num_adjacent, random_state=1
    )

    balanced_publications = pd.concat(
        [publications[publications["location"] == "Adjacent"], downstream_sample]
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


def _preprocess_data_fig2b(publications: pd.DataFrame) -> pd.DataFrame:

    # create a random sample of 10_000 publications of each source
    sample = (
        publications.groupby("source")
        .apply(lambda x: x.sample(10_000, random_state=42))
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
            "af": "AlphaFold",
            "ct_ai": "AI-intensive Frontier",
            "ct_noai": "Non-AI Frontier",
            "other": "Other Struct. Biol.",
        }
    )

    return field_counts


def _create_chart_fig2a(data: pd.DataFrame, title: str) -> alt.Chart:
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


def _create_chart_fig2b(data: pd.DataFrame, title: str) -> alt.Chart:
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


def generate_fig2(publications: pd.DataFrame) -> Image.Image:
    """Generate all figures 1a through 1d."""
    data_2a = _preprocess_data_fig2a(publications)
    data_2b = _preprocess_data_fig2b(publications)

    data_2a = data_2a[data_2a["primary_field"].isin(MAIN_TOPICS)]
    data_2b = data_2b[data_2b["primary_field"].isin(MAIN_TOPICS)]

    chart_a = _create_chart_fig2a(
        data_2a,
        "(a) Within-Field Shares by Distance to Core Research",
    )
    chart_b = _create_chart_fig2b(
        data_2b,
        "(b) Within-Field Shares by Core Source",
    )

    hconcat = alt.vconcat(chart_a, chart_b).resolve_scale(color="independent")

    image = save_chart_as_image(hconcat)

    return image


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
