"""
This is a boilerplate pipeline 'data_analysis_syllabii'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
import altair as alt

logger = logging.getLogger(__name__)

def compute_trends(
    data: pd.DataFrame,
):
    """
    Computes trends for syllabii based on the provided data.

    Args:
        data (pd.DataFrame): The input data containing syllabii information.

    Returns:
        tuple: A tuple containing two DataFrames:
            - trends: DataFrame with trends computed based on year and label.
            - region_trends: DataFrame with trends computed based on year, label, and region.
    """
    
    logger.info("Computing trends for syllabii")

    data.fillna({"label": "other", "country": "United States"}, inplace=True)
    # lets drop Taiwan, too few obs
    data = data[data["country"] != "Taiwan"]


    # Define country groups
    country_groups = {
        "United States": "United States",
        "United Kingdom": "United Kingdom",
        "Canada": "Canada",
        "Australia": "Australia & NZ",
        "New Zealand": "Australia & NZ",
        "Japan": "Japan"
    }

    # Replace countries with their groups
    data["region"] = data["country"].map(country_groups)
    data.fillna({"region": "Europe"}, inplace=True)


    trends = data.groupby(["year", "label"]).size().reset_index(name="count")
    region_trends = data.groupby(["year", "label", "region"]).size().reset_index(name="count")

    return trends, region_trends

def plot_trends(trends: pd.DataFrame, region_trends: pd.DataFrame):
    """
    Plots trends of the data.

    Args:
        trends (pd.DataFrame): The computed trends.
        region_trends (pd.DataFrame): The computed region trends.

    Returns:
        None
    """
    # Plot overall trends
    overall_chart = alt.Chart(trends).mark_line().encode(
        x="year:O",
        y="count:Q",
        color="label:N",
        tooltip=["year", "label", "count"]
    ).properties(
        title="Overall Trends per Label"
    )

    overall_chart.show()

    # Define regions
    regions = ["United States", "United Kingdom", "Canada", "Australia & NZ", "Japan", "Europe"]

    # Plot trends for each region
    region_charts = []
    for region in regions:
        region_data = region_trends[region_trends["region"] == region]
        chart = alt.Chart(region_data).mark_line().encode(
            x="year:O",
            y="count:Q",
            color="label:N",
            tooltip=["year", "label", "count"]
        ).properties(
            title=f"Trends per Label for {region}"
        )
        region_charts.append(chart)

    # Concatenate charts
    chart = alt.vconcat(*region_charts)

    chart.show()