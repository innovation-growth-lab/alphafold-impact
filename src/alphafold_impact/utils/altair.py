"""Altair utils"""

from io import BytesIO
import vl_convert as vlc
from PIL import Image
import altair as alt


def altair_to_png(altair_chart: alt.Chart) -> Image.Image:
    """Convert Altair chart to PNG file"""
    png_data = vlc.vegalite_to_png(altair_chart.to_json(), scale=2)
    return Image.open(BytesIO(png_data))
