"""
This is a boilerplate pipeline 'data_processing_os'
generated using Kedro 0.19.1
"""

import logging
import pandas as pd
from ...utils import scrape_constants  # pylint: disable=relative-beyond-top-level

EXPERIMENTAL_KEYWORDS = scrape_constants.EXPERIMENTAL_KEYWORDS
EXPERIMENTAL_LOOSE_KEYWORDS = scrape_constants.EXPERIMENTAL_LOOSE_KEYWORDS
COMPUTATIONAL_KEYWORDS = scrape_constants.COMPUTATIONAL_KEYWORDS
STRUCTURAL_BIOLOGY_KEYWORDS = scrape_constants.STRUCTURAL_BIOLOGY_KEYWORDS

logger = logging.getLogger(__name__)


def get_regional_biology_syllabii(data_loaders):
    """
    Retrieves regional biology syllabii from data loaders.

    Args:
        data_loaders (dict): A dictionary of data loaders.

    Returns:
        pandas.DataFrame: A concatenated DataFrame containing the regional biology syllabii.
    """
    outputs = []
    for i, data_loader in enumerate(data_loaders.values()):
        logger.info("Loading data %d / %d", i + 1, len(data_loaders))
        outputs.append(data_loader())
    outputs = pd.concat(outputs, ignore_index=True)
    outputs["country"].fillna("United States", inplace=True)
    outputs["country"] = outputs["country"].str.replace("+", " ")
    return outputs


def get_topic_biology_syllabii(data_loaders):
    """
    Retrieves topic, label, and keyword information from data loaders.

    Args:
        data_loaders (dict): A dictionary of data loaders.

    Returns:
        pandas.DataFrame: A DataFrame containing the retrieved information.
    """
    outputs = []
    for i, (keyword, data_loader) in enumerate(data_loaders.items()):
        logger.info("Loading data %d / %d", i + 1, len(data_loaders))
        data = data_loader()
        keyword = keyword[9:]
        if "_" in keyword:
            keyword = keyword.split("_")[0]
        data["keyword"] = keyword

        # add label for experimental, computational, or structural biology
        if keyword in EXPERIMENTAL_KEYWORDS:
            data["label"] = "experimental"
        elif keyword in EXPERIMENTAL_LOOSE_KEYWORDS:
            data["label"] = "experimental (loose)"
        elif keyword in COMPUTATIONAL_KEYWORDS:
            data["label"] = "computational"
        elif keyword in STRUCTURAL_BIOLOGY_KEYWORDS:
            data["label"] = "structural biology"
        else:
            data["label"] = "other"

        outputs.append(data)

    outputs = pd.concat(outputs, ignore_index=True)

    return outputs


def label_biology_courses(regional_data, keyword_data):
    """
    Labels biology courses based on regional and keyword data.

    Args:
        regional_data (pandas.DataFrame): The regional data containing course
            information.
        keyword_data (pandas.DataFrame): The keyword data containing course 
            labels and text snippets.

    Returns:
        pandas.DataFrame: The merged and labeled data.

    """
    # merge regional and keyword data
    merged_data = pd.merge(
        regional_data,
        keyword_data[["course_code", "label", "text_snippet"]],
        on="course_code",
        how="left",
    )

    # drop duplicate rows based on course_code, year, label
    merged_data.drop_duplicates(subset=["course_code", "year", "label"], inplace=True)

    # reset index
    merged_data.reset_index(drop=True, inplace=True)

    merged_data["label"].fillna("other", inplace=True)

    return merged_data
