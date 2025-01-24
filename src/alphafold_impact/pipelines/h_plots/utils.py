import logging
from io import BytesIO
import pandas as pd
import altair as alt
import vl_convert as vlc
from PIL import Image

COUNTRY_CLASSIFICATION = {
    "AE": "Non-LMIC",
    "AG": "Non-LMIC",
    "AL": "LMIC",
    "AM": "LMIC",
    "AO": "LMIC",
    "AR": "Non-LMIC",
    "AT": "Non-LMIC",
    "AU": "Non-LMIC",
    "AW": "Non-LMIC",
    "AZ": "LMIC",
    "BA": "LMIC",
    "BB": "Non-LMIC",
    "BD": "LMIC",
    "BE": "Non-LMIC",
    "BF": "LMIC",
    "BG": "LMIC",
    "BH": "Non-LMIC",
    "BI": "LMIC",
    "BJ": "LMIC",
    "BN": "Non-LMIC",
    "BO": "LMIC",
    "BR": "Non-LMIC",
    "BS": "Non-LMIC",
    "BT": "LMIC",
    "BW": "LMIC",
    "BY": "LMIC",
    "BZ": "Non-LMIC",
    "CA": "Non-LMIC",
    "CD": "LMIC",
    "CF": "LMIC",
    "CG": "Non-LMIC",
    "CH": "Non-LMIC",
    "CI": "Non-LMIC",
    "CL": "Non-LMIC",
    "CM": "LMIC",
    "CN": "LMIC",
    "CO": "LMIC",
    "CR": "LMIC",
    "CU": "Non-LMIC",
    "CY": "Non-LMIC",
    "CZ": "Non-LMIC",
    "DE": "Non-LMIC",
    "DK": "Non-LMIC",
    "DM": "LMIC",
    "DO": "Non-LMIC",
    "DZ": "LMIC",
    "EC": "LMIC",
    "EE": "Non-LMIC",
    "EG": "LMIC",
    "ES": "Non-LMIC",
    "ET": "LMIC",
    "FI": "Non-LMIC",
    "FJ": "LMIC",
    "FO": "Non-LMIC",
    "FR": "Non-LMIC",
    "GB": "Non-LMIC",
    "GD": "LMIC",
    "GE": "LMIC",
    "GF": "Non-LMIC",
    "GH": "LMIC",
    "GM": "LMIC",
    "GP": "Non-LMIC",
    "GR": "Non-LMIC",
    "GT": "LMIC",
    "GU": "Non-LMIC",
    "HK": "Non-LMIC",
    "HR": "Non-LMIC",
    "HU": "Non-LMIC",
    "ID": "LMIC",
    "IE": "Non-LMIC",
    "IL": "Non-LMIC",
    "IN": "LMIC",
    "IQ": "LMIC",
    "IR": "LMIC",
    "IS": "Non-LMIC",
    "IT": "Non-LMIC",
    "JM": "LMIC",
    "JO": "LMIC",
    "JP": "Non-LMIC",
    "KE": "LMIC",
    "KG": "LMIC",
    "KH": "LMIC",
    "KN": "Non-LMIC",
    "KR": "Non-LMIC",
    "KW": "Non-LMIC",
    "KZ": "LMIC",
    "LB": "LMIC",
    "LK": "LMIC",
    "LS": "LMIC",
    "LT": "Non-LMIC",
    "LU": "Non-LMIC",
    "LV": "Non-LMIC",
    "LY": "Non-LMIC",
    "MA": "LMIC",
    "MC": "Non-LMIC",
    "MD": "Non-LMIC",
    "ME": "LMIC",
    "MG": "LMIC",
    "MK": "LMIC",
    "ML": "LMIC",
    "MM": "LMIC",
    "MN": "LMIC",
    "MO": "Non-LMIC",
    "MT": "Non-LMIC",
    "MU": "Non-LMIC",
    "MW": "LMIC",
    "MX": "Non-LMIC",
    "MY": "Non-LMIC",
    "MZ": "LMIC",
    "NC": "Non-LMIC",
    "NE": "LMIC",
    "NG": "LMIC",
    "NI": "LMIC",
    "NL": "Non-LMIC",
    "NO": "Non-LMIC",
    "NP": "LMIC",
    "NZ": "Non-LMIC",
    "OM": "Non-LMIC",
    "PA": "Non-LMIC",
    "PE": "LMIC",
    "PF": "Non-LMIC",
    "PG": "LMIC",
    "PH": "LMIC",
    "PK": "LMIC",
    "PL": "Non-LMIC",
    "PR": "Non-LMIC",
    "PS": "LMIC",
    "PT": "Non-LMIC",
    "PY": "Non-LMIC",
    "QA": "Non-LMIC",
    "RE": "Non-LMIC",
    "RO": "Non-LMIC",
    "RS": "LMIC",
    "RU": "Non-LMIC",
    "RW": "LMIC",
    "SA": "Non-LMIC",
    "SD": "LMIC",
    "SE": "Non-LMIC",
    "SG": "Non-LMIC",
    "SI": "Non-LMIC",
    "SK": "Non-LMIC",
    "SL": "LMIC",
    "SN": "LMIC",
    "SS": "LMIC",
    "SY": "LMIC",
    "TG": "LMIC",
    "TH": "LMIC",
    "TJ": "LMIC",
    "TN": "LMIC",
    "TR": "Non-LMIC",
    "TT": "LMIC",
    "TW": "Non-LMIC",
    "TZ": "LMIC",
    "UA": "LMIC",
    "UG": "LMIC",
    "US": "Non-LMIC",
    "UY": "Non-LMIC",
    "UZ": "LMIC",
    "VE": "LMIC",
    "VI": "Non-LMIC",
    "VN": "LMIC",
    "XK": "LMIC",
    "YE": "LMIC",
    "ZA": "LMIC",
    "ZM": "LMIC",
    "ZW": "LMIC",
}

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
