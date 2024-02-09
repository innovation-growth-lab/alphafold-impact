"""
This module contains functions for processing NIH funding opportunities data.

Functions:
    clean_nih_funding_opportunites: Cleans NIH funding opportunities data.
"""

import pandas as pd


def _remove_leading_text_from_cols(
    df: pd.DataFrame, leading_text_to_remove: str
) -> pd.DataFrame:
    """
    Removes specified leading text from column names in a pandas DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame whose columns need to be renamed.
        leading_text_to_remove (str): The leading text to be removed from the column names.

    Returns:
        pd.DataFrame: A DataFrame with specified leading text removed from column names.

    """
    return df.rename(
        columns=lambda col: (
            col[len(leading_text_to_remove) :]
            if col.startswith(leading_text_to_remove)
            else col
        )
    )


def _filter_dataframe_by_year(
    df: pd.DataFrame, column_name: str, year: int
) -> pd.DataFrame:
    """
    Filter a DataFrame based on a year condition in a specified datetime column.

    Args:
        df (pd.DataFrame): The DataFrame to be filtered.
        column_name (str): The name of the column containing datetime information.
        year (int): The year to include records from.

    Returns:
        pd.DataFrame: A DataFrame filtered to include rows from the specified year onwards.

    Example:
        >>> df = pd.DataFrame({'reldate': ['1997-02-07T00:00:00.000Z', '1995-01-01T00:00:00.000Z']})
        >>> filter_dataframe_by_year(df, 'reldate')
        # returns df with only the row '1997-02-07T00:00:00.000Z'
    """
    df[column_name] = pd.to_datetime(df[column_name])
    return df[df[column_name].dt.year >= year]


def clean_nih_funding_opportunites(
    nih_funding_opportunities_raw: pd.DataFrame, year: int
) -> pd.DataFrame:
    """
    Clean NIH funding opportunities data by updating column names
    and filtering the data to include records from the specified year

    Args:
        nih_funding_opportunities_raw (pd.DataFrame): NIH funding opportunities data
        year (int): The year to include records from.


    Returns:
        pd.DataFrame: Cleaned NIH funding opportunities data
    """
    return (
        nih_funding_opportunities_raw.pipe(_remove_leading_text_from_cols, "_")
        .pipe(_remove_leading_text_from_cols, "source.")
        .pipe(_filter_dataframe_by_year, "reldate", year)
        .drop(columns=["index"])
    )
