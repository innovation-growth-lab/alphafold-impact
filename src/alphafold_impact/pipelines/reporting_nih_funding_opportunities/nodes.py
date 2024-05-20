"""This module contains functions to generate charts for NIH funding opportunities."""

import pandas as pd
import altair as alt
from PIL import Image
from alphafold_impact.utils.altair import altair_to_png


def _prepare_fadirectcosts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts the 'fadirectcosts' column to numeric and drops rows with NaN values in that column.

    Args:
        df (pd.DataFrame): The DataFrame containing a 'fadirectcosts' column.

    Returns:
        pd.DataFrame: A DataFrame with 'fadirectcosts' as numeric and NaNs removed.
    """
    df["fadirectcosts"] = pd.to_numeric(df["fadirectcosts"], errors="coerce")
    return df.dropna(subset=["fadirectcosts"])


def plot_eb_sb_funding_opportunities(
    nih_sb_funding_opps: pd.DataFrame, nih_eb_funding_opps: pd.DataFrame
) -> Image.Image:
    """
    Generate a stacked bar chart of counts of funding opportunities per year for the
    different subfields of Structural Biology and Experimental Biology.

    Args:
        nih_sb_funding_opps (pd.DataFrame): NIH Structural Biology funding opportunities DataFrame
            with 'reldate', and 'subfield' columns.
        nih_eb_funding_opps (pd.DataFrame): NIH Experimental Biology funding opportunities DataFrame
            with 'reldate', and 'subfield' columns.

    Returns:
        Image.Image: PNG chart.
    """
    combined_df = pd.concat([nih_sb_funding_opps, nih_eb_funding_opps])
    combined_df["year"] = pd.to_datetime(combined_df["reldate"]).dt.year
    grouped = combined_df.groupby(["year", "subfield"]).size().reset_index(name="count")
    eb_sb_funding_opportunities = (
        alt.Chart(grouped)
        .mark_bar()
        .encode(
            x=alt.X("year:O", axis=alt.Axis(title="Year", labelAngle=0)),
            y=alt.Y(
                "count:Q",
                stack="zero",
                axis=alt.Axis(title="Funding Opportunities"),
            ),
            color=alt.Color(
                "subfield:N",
                scale=alt.Scale(range=["#1f77b4", "#aec7e8"]),
                legend=alt.Legend(
                    orient="none",
                    title="Subfield",
                    legendX=10,
                    legendY=10,
                    padding=10,
                ),
            ),
            tooltip=["year", "subfield", "count"],
        )
        .properties(title="SB and EB NIH Funding Opportunities per Year", width=600)
    )
    return altair_to_png(eb_sb_funding_opportunities)


def plot_eb_sb_funding_opportunity_amounts(
    nih_sb_funding_opps: pd.DataFrame, nih_eb_funding_opps: pd.DataFrame
) -> Image.Image:
    """
    Generate a stacked bar chart of counts of funding opportunity amounts per year for the
    different subfields of Structural Biology and Experimental Biology.

    Args:
        nih_sb_funding_opps (pd.DataFrame): NIH Structural Biology funding opportunities DataFrame
            with 'reldate', 'subfield' and 'fadirectcosts' columns.
        nih_eb_funding_opps (pd.DataFrame): NIH Experimental Biology funding opportunities DataFrame
            with 'reldate', 'subfield' and 'fadirectcosts' columns.

    Returns:
        Image.Image: PNG chart.
    """
    combined_df = pd.concat(
        [
            _prepare_fadirectcosts(nih_sb_funding_opps),
            _prepare_fadirectcosts(nih_eb_funding_opps),
        ]
    )
    combined_df["year"] = pd.to_datetime(combined_df["reldate"]).dt.year
    grouped = (
        combined_df.groupby(["year", "subfield"])["fadirectcosts"].sum().reset_index()
    )
    eb_sb_funding_opportunity_amounts = (
        alt.Chart(grouped)
        .mark_bar()
        .encode(
            x=alt.X("year:O", axis=alt.Axis(title="Year", labelAngle=0)),
            y=alt.Y(
                "fadirectcosts:Q",
                stack="zero",
                axis=alt.Axis(title="Funding Opportunity Amounts ($)"),
            ),
            color=alt.Color(
                "subfield:N",
                scale=alt.Scale(range=["#1f77b4", "#aec7e8"]),
                legend=alt.Legend(
                    orient="none",
                    title="Subfield",
                    legendX=10,
                    legendY=10,
                    padding=10,
                ),
            ),
            tooltip=["year", "subfield", "fadirectcosts"],
        )
        .properties(
            title="SB and EB NIH Funding Opportunity Amounts per Year", width=600
        )
    )
    return altair_to_png(eb_sb_funding_opportunity_amounts)


def plot_total_funding_opportunities(
    nih_funding_opps: pd.DataFrame,
) -> Image.Image:
    """
    Generate a bar chart of counts of total funding opportunities per year.

    Args:
        nih_funding_opps (pd.DataFrame): NIH funding opportunities DataFrame
            containing funding opportunities with 'reldate' column.

    Returns:
        Image.Image: PNG chart.
    """
    nih_funding_opps["year"] = pd.to_datetime(nih_funding_opps["reldate"]).dt.year
    grouped = nih_funding_opps.groupby(["year"]).size().reset_index(name="count")
    total_funding_opps = (
        alt.Chart(grouped)
        .mark_bar()
        .encode(
            x=alt.X("year:O", axis=alt.Axis(title="Year", labelAngle=0)),
            y=alt.Y(
                "count:Q",
                stack="zero",
                axis=alt.Axis(title="Funding Opportunities"),
            ),
            tooltip=["year", "count"],
        )
        .properties(title="Total NIH Funding Opportunities", width=600)
    )
    return altair_to_png(total_funding_opps)


def plot_total_funding_opportunity_amounts(
    nih_funding_opps: pd.DataFrame,
) -> Image.Image:
    """
    Generate a bar chart of counts of total funding opportunity amounts per
    year.

    Args:
        nih_funding_opps (pd.DataFrame): NIH funding oppportunities DataFrame
            which contains funding opportunities with 'reldate' and
            'fadirectcosts' columns.

    Returns:
        Image.Image: PNG chart.
    """
    nih_funding_opps["year"] = pd.to_datetime(nih_funding_opps["reldate"]).dt.year
    grouped = (
        _prepare_fadirectcosts(nih_funding_opps)
        .groupby(["year"])["fadirectcosts"]
        .sum()
        .reset_index()
    )
    total_funding_opp_amounts_per_year = (
        alt.Chart(grouped)
        .mark_bar()
        .encode(
            x=alt.X("year:O", axis=alt.Axis(title="Year", labelAngle=0)),
            y=alt.Y(
                "fadirectcosts:Q",
                stack="zero",
                axis=alt.Axis(title="Funding Opportunity Amounts ($)"),
            ),
            tooltip=["year", "fadirectcosts"],
        )
        .properties(title="Total NIH Funding Opportunity Amounts", width=600)
    )
    return altair_to_png(total_funding_opp_amounts_per_year)
