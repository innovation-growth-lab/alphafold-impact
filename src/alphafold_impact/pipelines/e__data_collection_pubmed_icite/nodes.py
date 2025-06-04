"""
This module contains functions for fetching ClinicalTrials.gov data using API v2.

Functions:
    collect_and_save_clinical_trials_v2: Downloads, extracts, combines
        and saves ClinicalTrials.gov data as a partitioned dataset using the new API v2.
"""

import logging
from pathlib import Path
from typing import List, Dict, Callable
import zipfile
import shutil
import contextlib
import json

import requests
import pandas as pd

logger = logging.getLogger(__name__)


def _download_zip_file(zip_url: str, zip_save_path: str) -> Path:
    """
    Download a ZIP file from the new API v2 endpoint.

    Args:
        zip_url (str): URL of the zip file to be downloaded.
        zip_save_path (str): Path to save the downloaded file.

    Returns:
        Path: The path to the downloaded zip file.
    """
    logger.info("Downloading ZIP file from API v2 endpoint: %s", zip_url)

    # Add headers for the API request
    headers = {
        "User-Agent": "AlphaFold-Impact-Research/1.0",
        "Accept": "application/zip",
    }

    response = requests.get(
        zip_url, headers=headers, timeout=300
    )  # Increased timeout for large files
    response.raise_for_status()

    zip_path = Path(zip_save_path)
    zip_path.parent.mkdir(parents=True, exist_ok=True)

    with open(zip_path, "wb") as file:
        file.write(response.content)

    logger.info(
        "Downloaded ZIP file to %s (size: %d bytes)", zip_path, zip_path.stat().st_size
    )
    return zip_path


def _unzip_and_flatten(zip_path: Path, extract_path: str) -> List[Path]:
    """
    Unzip a file and move all contents to a main directory, flattening the structure.
    The new API v2 provides JSON files instead of XML.

    Args:
        zip_path (Path): Path to the zip file.
        extract_path (str): Path to extract and store the contents of the zip file.

    Returns:
        List[Path]: List of paths to the flattened JSON data files.
    """
    extract_path = Path(extract_path)
    logger.info("Unzipping API v2 ZIP file %s to %s", zip_path, extract_path)
    flattened_file_paths = []

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(extract_path)

    logger.info("Flattening folder structure and collecting JSON files")
    for path in extract_path.rglob("*.json"):  # Only process JSON files
        if path.is_file() and path.name != "Contents.txt":  # Skip the contents file
            destination = extract_path / path.name
            if path != destination:
                shutil.move(str(path), str(destination))
            flattened_file_paths.append(destination)

    logger.info("Deleting empty folders")
    for path in extract_path.iterdir():
        if path.is_dir():
            with contextlib.suppress(OSError):
                shutil.rmtree(path)  # Use rmtree for non-empty directories

    logger.info(
        "Unzipped and flattened %d JSON files to %s",
        len(flattened_file_paths),
        extract_path,
    )

    zip_path.unlink()
    logger.info("Deleted ZIP file %s", zip_path)

    return flattened_file_paths


def _load_and_normalize_json_v2(file_path: Path) -> pd.DataFrame:
    """
    Loads ClinicalTrials.gov API v2 JSON file and normalizes the study data
    into a pandas DataFrame. The new API v2 has a simplified structure.

    Args:
        file_path (Path): The file path of the JSON file.

    Returns:
        pd.DataFrame: A DataFrame with the normalized data.
    """
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            json_data = json.load(file)

        # The new API v2 structure is much flatter
        # Each JSON file contains a single study record
        if isinstance(json_data, dict):
            # Normalize the study data - the structure is already flatter in v2
            df = pd.json_normalize(json_data)
            return df
        else:
            logger.warning("Unexpected JSON structure in file %s", file_path)
            return pd.DataFrame()

    except (json.JSONDecodeError, KeyError) as e:
        logger.error("Error processing file %s: %s", file_path, e)
        return pd.DataFrame()


def _process_batch_v2(files_to_process: List[Path], batch_num: str) -> pd.DataFrame:
    """
    Processes a batch of JSON files from API v2, concatenating their contents into
    a single DataFrame.

    This function loads JSON data from each file in the batch,
    normalizes it, and combines it into a single DataFrame.
    After processing, each JSON file is deleted.

    Args:
        files_to_process (List[Path]): A list of file paths to process.
        batch_num (str): The batch number, used for logging purposes.

    Returns:
        pd.DataFrame: A DataFrame containing the combined data from all files in the batch.
    """
    logger.info(
        "Processing API v2 batch %s with %d files", batch_num, len(files_to_process)
    )

    dataframes = []
    for fp in files_to_process:
        df = _load_and_normalize_json_v2(fp)
        if not df.empty:
            dataframes.append(df)

    if dataframes:
        combined_df = pd.concat(dataframes, sort=False, ignore_index=True)
        logger.info("Combined %d studies in batch %s", len(combined_df), batch_num)
    else:
        combined_df = pd.DataFrame()
        logger.warning("No valid data found in batch %s", batch_num)

    # Delete files after processing
    for file_path in files_to_process:
        file_path.unlink()

    logger.info("Deleted JSON files related to batch %s", batch_num)
    logger.info("Completed processing batch %s", batch_num)
    return combined_df


def _json_to_csv_v2(
    dir_path: str, file_paths: List[Path], batch_size: int = 1000
) -> Dict[str, Callable]:
    """
    Processes JSON files from API v2 in the given directory and saves them in batches
    of specified records as CSV files. Each batch processing is deferred
    until the callable is executed.

    Args:
        dir_path (str): Path to the directory containing JSON files.
        file_paths (List[Path]): List of JSON file paths to be processed.
        batch_size (int): Number of files (studies) in each batch.
                         Reduced from 20000 since each file is now one study.

    Returns:
        Dict[str, Callable]: A dictionary where each key is the name of
        a CSV file, and each value is a callable that generates and
        returns the DataFrame to be saved in that CSV file.
    """
    logger.info("Starting API v2 JSON to CSV processing for %d files", len(file_paths))
    dir_path = Path(dir_path)
    partitions = {}

    for batch_number, i in enumerate(range(0, len(file_paths), batch_size), start=1):
        batch_files = file_paths[i : i + batch_size]
        current_batch_number = str(batch_number).zfill(3)
        file_name = f"clinical_trials_v2_batch_{current_batch_number}"
        partitions[file_name] = (
            lambda files_to_process=batch_files, batch_num=current_batch_number: _process_batch_v2(
                files_to_process, batch_num
            )
        )

    logger.info(
        "Completed setting up callables for API v2 JSON to CSV processing (%d batches)",
        len(partitions),
    )
    return partitions


def collect_and_save_clinical_trials_v2(
    api_v2_download_url: str, zip_save_path: str, extract_path: str
) -> Dict[str, Callable]:
    """
    Downloads, extracts, combines and saves ClinicalTrials.gov data using API v2.

    Args:
        api_v2_download_url (str): URL of the API v2 studies download endpoint
            (should be /api/v2/studies/download?format=json.zip).
        zip_save_path (str): Path to save the downloaded ZIP file.
        extract_path (str): Path to extract and store the
            contents of the ZIP file.

    Returns:
        Dict[str, Callable]: A dictionary with keys being the CSV file
        names and values being callables that return DataFrames
        representing the data for each CSV file.
    """
    logger.info("Starting ClinicalTrials.gov API v2 data collection")

    zip_path = _download_zip_file(api_v2_download_url, zip_save_path)
    json_file_paths = _unzip_and_flatten(zip_path, extract_path)

    if not json_file_paths:
        logger.warning("No JSON files found in the downloaded archive")
        return {}

    return _json_to_csv_v2(extract_path, json_file_paths)


# Keep the legacy function for backward compatibility during transition
def collect_and_save_nih_clinical_trials(
    nih_ct_download_url: str, nih_ct_zip_save_path: str, nih_ct_extract_path: str
) -> Dict[str, Callable]:
    """
    Legacy function - redirects to the new API v2 implementation.

    This function is kept for backward compatibility during the migration period.
    It will call the new API v2 function instead of the old implementation.

    Args:
        nih_ct_download_url (str): URL of the clinical trials data download.
        nih_ct_zip_save_path (str): Path to save the downloaded ZIP file.
        nih_ct_extract_path (str): Path to extract and store the contents.

    Returns:
        Dict[str, Callable]: A dictionary with keys being the CSV file
        names and values being callables that return DataFrames.
    """
    logger.warning(
        "Using legacy function name - consider updating to 'collect_and_save_clinical_trials_v2'"
    )
    return collect_and_save_clinical_trials_v2(
        nih_ct_download_url, nih_ct_zip_save_path, nih_ct_extract_path
    )
