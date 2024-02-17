"""
This file contains functions for using the SKR web python API
for submitting batch jobs using the NLM Medical Text Indexer (MTI).
It can be used to add MeSH labels to text.

See here for installation steps: https://github.com/lhncbc/skr_web_python_api

Get UMLS API key from here: https://uts.nlm.nih.gov/uts/profile
"""

from skr_web_api import Submission
import pandas as pd


def _remove_newlines(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """
    Remove newline characters from a specified column in a pandas DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        column (str): The name of the column from which to remove newlines.

    Returns:
        pd.DataFrame: The DataFrame with newline characters removed from the specified column.
    """
    df[column] = df[column].astype(str).str.replace("\n", "", regex=False)
    return df


def _remove_prc_rows(mesh_results: str) -> str:
    """
    Processes the input string, removing lines that contain '|PRC|'.

    Args:
        mesh_results (str): A string containing the input data.

    Returns:
        str: A string with lines containing '|PRC|' removed.
    """
    return "\n".join(
        [line for line in mesh_results.splitlines() if "|PRC|" not in line]
    )


def generate_text_for_mesh_tagging(df: pd.DataFrame, id_col: str, text_col: str) -> str:
    """
    Generates a formatted string into a suitable format for MeSH tagging.

    Args:
        df (pd.DataFrame): The DataFrame to process.
        id_col (str): The name of the column in the DataFrame to use for IDs.
        text_col (str): The name of the column in the DataFrame to use for text.

    Returns:
        str: A string containing the formatted text, with each line in the format f'{id}|{text}\n'.
    """
    df_text_newlines_removed = _remove_newlines(df, text_col)
    lines = [
        f"{row[id_col]}|{row[text_col]}"
        for _, row in df_text_newlines_removed.iterrows()
    ]
    return "\n".join(lines)


def skr_web_python_api_generic_batch(
    email: str,
    api_key: str,
    batch_name: str,
    text_to_tag: str,
    cmd: str,
    cmdargs: str,
) -> str:
    """Submits NIH Generic Batch job using the SKR Web API

    Args:
        email (str): User's email address.
        api_key (str): UTS API key.
        batch_name (str): Name to give to batch job.
        text_to_tag (str): Text to submit to be tagged. Each record
            should be on a new line and be in the format {id}|{text}.
            See example here: https://lhncbc.nlm.nih.gov/ii/about/help/single_line_ID_sample.txt
        cmd (str): Command to be executed.
        cmdargs (str): Arguments for the command. A full list of args can be found:
            https://lhncbc.nlm.nih.gov/ii/tools/MTI/help_info.html
    """
    # Move this to be a parameter ####
    serviceurl = (
        "https://ii.nlm.nih.gov/cgi-bin/II/UTS_Required/API_batchValidationII.pl"
    )

    inst = Submission(email, api_key)
    inst.set_serviceurl(serviceurl)
    inst.init_generic_batch(cmd, cmdargs)
    inst.set_batch_file(batch_name, text_to_tag)
    inst.form["SingLinePMID"] = True
    response = inst.submit()
    return response.content.decode().replace("NOT DONE LOOP\n", "")


def mesh_results_to_df(mesh_results: str) -> pd.DataFrame:
    """
    Converts the input string data into a DataFrame with
    columns 'text_id', 'mesh_term', 'mesh_term_id'.

    Args:
        input_data (str): A string containing the input data.

    Returns:
        pd.DataFrame: A DataFrame with columns, 'text_id',
            'mesh_term', 'mesh_term_id'.
    """
    mesh_results_no_prc = _remove_prc_rows(mesh_results)
    rows = [line.split("|") for line in mesh_results_no_prc.strip().split("\n")]
    return pd.DataFrame(
        rows, columns=["text_id", "mesh_term", "mesh_term_id", "position"]
    ).drop(columns=["position"])
