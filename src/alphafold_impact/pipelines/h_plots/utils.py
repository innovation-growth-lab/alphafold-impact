import logging
from io import BytesIO
import pandas as pd
import altair as alt
import vl_convert as vlc
from PIL import Image

logger = logging.getLogger(__name__)


def save_chart_as_image(chart: alt.Chart) -> Image.Image:
    """Convert an Altair chart to a Pillow Image object."""
    png_bytes = vlc.vegalite_to_png(  # pylint: disable=no-member
        chart.to_json(), scale=3
    )
    image = Image.open(BytesIO(png_bytes))
    return image


def normalise_january_counts(monthly_counts: pd.DataFrame):
    """Normalise January counts in the input DataFrame."""
    adjusted_counts = monthly_counts.copy()

    adjusted_counts["count"] = adjusted_counts["count"].astype(float)

    for (year, adjacent, downstream), group in monthly_counts.groupby(
        ["year", "adjacent", "downstream"]
    ):
        if 1 in group["month"].values:
            december_count = monthly_counts[
                (monthly_counts["year"] == year - 1)
                & (monthly_counts["month"] == 12)
                & (monthly_counts["adjacent"] == adjacent)
                & (monthly_counts["downstream"] == downstream)
            ]["count"].sum()
            february_count = monthly_counts[
                (monthly_counts["year"] == year)
                & (monthly_counts["month"] == 2)
                & (monthly_counts["adjacent"] == adjacent)
                & (monthly_counts["downstream"] == downstream)
            ]["count"].sum()

            # impute January count as mean of December and February
            imputed_january_count = (december_count + february_count) / 2

            # calculate original January count
            original_january_count = group[group["month"] == 1]["count"].values[0]

            # compute the adjustment difference
            adjustment_difference = original_january_count - imputed_january_count

            # calculate the number of months in the current group
            num_months = group["month"].nunique()

            # apply the adjustment across all months
            monthly_adjustment = adjustment_difference / num_months

            # update January count
            adjusted_counts.loc[
                (adjusted_counts["year"] == year)
                & (adjusted_counts["month"] == 1)
                & (adjusted_counts["adjacent"] == adjacent)
                & (adjusted_counts["downstream"] == downstream),
                "count",
            ] = imputed_january_count

            # adjust all months
            adjusted_counts.loc[
                (adjusted_counts["year"] == year)
                & (adjusted_counts["adjacent"] == adjacent)
                & (adjusted_counts["downstream"] == downstream),
                "count",
            ] += monthly_adjustment

    return adjusted_counts
