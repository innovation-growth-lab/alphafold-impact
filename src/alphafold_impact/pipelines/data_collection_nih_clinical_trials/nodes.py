"""
This module contains functions for fetching NIH clinical trials data.

Functions:
    collect_and_save_nih_clinical_trials: Downloads, extracts, combines
        and saves NIH clinical trials data as a partitioned dataset.
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
    Download a ZIP file from a given URL.

    Args:
        zip_url (str): URL of the zip file to be downloaded.
        zip_save_path (str): Path to save the downloaded file.

    Returns:
        Path: The path to the downloaded zip file.
    """
    logger.info("Downloading ZIP file from %s", zip_url)
    response = requests.get(zip_url, timeout=30)
    response.raise_for_status()
    zip_save_path = Path(zip_save_path)
    zip_save_path.parent.mkdir(parents=True, exist_ok=True)
    with open(zip_save_path, "wb") as file:
        file.write(response.content)

    logger.info("Downloaded ZIP file to %s", zip_save_path)
    return zip_save_path


def _unzip_and_flatten(zip_path: Path, extract_path: str) -> List[Path]:
    """
    Unzip a file and move all contents to a main directory, flattening the structure.

    Args:
        zip_path (Path): Path to the zip file.
        extract_path (str): Path to extract and store the contents of the zip file.

    Returns:
        List[Path]: List of paths to the flattened data.
    """
    extract_path = Path(extract_path)
    logger.info("Unzipping ZIP file %s to %s", zip_path, extract_path)
    flattened_file_paths = []
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(extract_path)

    logger.info("Flattening folder structure")
    for path in extract_path.rglob("*"):
        if path.is_file():
            destination = extract_path / path.name
            if path != destination:
                shutil.move(str(path), str(destination))
                flattened_file_paths.append(destination)

    logger.info("Deleting empty folders")
    for path in extract_path.iterdir():
        if path.is_dir():
            with contextlib.suppress(OSError):
                path.rmdir()

    logger.info("Unzipped and flattened files to %s", extract_path)

    zip_path.unlink()

    logger.info("Deleted ZIP file %s", zip_path)

    return flattened_file_paths


def _load_and_normalize_json(file_path: Path) -> pd.DataFrame:
    """
    Loads NIH clinical trial JSON file and normalizes the
    'Study' data into a pandas DataFrame. The 'Rank' data
    is also added.

    Args:
        file_path (Path): The file path of the JSON file.

    Returns:
        pd.DataFrame: A DataFrame with the normalized data.
    """
    with open(file_path, "r", encoding="utf-8") as file:
        json_data = json.load(file)

    return pd.json_normalize(json_data["FullStudy"]["Study"]).assign(
        Rank=json_data["FullStudy"]["Rank"]
    )


def _process_batch(files_to_process: List[Path], batch_num: str) -> pd.DataFrame:
    """
    Processes a batch of JSON files, concatenating their contents into
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
    logger.info("Processing batch %s with %d files", batch_num, len(files_to_process))
    dataframes = [_load_and_normalize_json(fp) for fp in files_to_process]
    combined_df = pd.concat(dataframes, sort=False, join="outer")

    # Delete files after processing
    for file_path in files_to_process:
        file_path.unlink()

    logger.info("Deleted JSON files related to batch %s", batch_num)
    logger.info("Completed processing batch %s", batch_num)
    return combined_df


def _json_to_csv(
    dir_path: str, file_paths: List[Path], batch_size: int = 20000
) -> Dict[str, Callable]:
    """
    Processes JSON files in the given directory and saves them in batches
    of specified records as CSV files. Each batch processing is deferred
    until the callable is executed.

    Args:
        dir_path (str): Path to the directory containing JSON files.
        file_paths (List[Path]): List of JSON file paths to be processed.
        batch_size (int): Number of records in each batch.

    Returns:
        Dict[str, Callable]: A dictionary where each key is the name of
        a CSV file, and each value is a callable that generates and
        returns the DataFrame to be saved in that CSV file.
    """
    logger.info("Starting JSON to CSV processing")
    dir_path = Path(dir_path)
    partitions = {}
    for batch_number, i in enumerate(range(0, len(file_paths), batch_size), start=1):
        batch_files = file_paths[i : i + batch_size]
        current_batch_number = str(batch_number).zfill(3)
        file_name = f"nih_clinical_trials_batch_{current_batch_number}.csv"
        partitions[
            file_name
        ] = lambda files_to_process=batch_files, batch_num=current_batch_number: _process_batch(
            files_to_process, batch_num
        )
    logger.info("Completed setting up callables for JSON to CSV processing")
    return partitions


def collect_and_save_nih_clinical_trials(
    nih_ct_download_url: str, nih_ct_zip_save_path: str, nih_ct_extract_path: str
) -> Dict[str, Callable]:
    """
    Downloads, extracts, combines and saves NIH clinical trials data.

    Args:
        nih_ct_download_url (str): URL of the NIH clinical trials ZIP
            file to be downloaded.
        nih_ct_zip_save_path (str): Path to save the downloaded ZIP file.
        nih_ct_extract_path (str): Path to extract and store the
            contents of the ZIP file.

    Returns:
        Dict[str, Callable]: A dictionary with keys being the CSV file
        names and values being callables that return DataFrames
        representing the data for each CSV file.
    """
    zip_path = _download_zip_file(nih_ct_download_url, nih_ct_zip_save_path)
    json_file_paths = _unzip_and_flatten(zip_path, nih_ct_extract_path)
    return _json_to_csv(nih_ct_extract_path, json_file_paths)
