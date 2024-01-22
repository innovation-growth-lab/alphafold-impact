# Gateway to Research (GtR) Data Collection Pipeline (`data_collection_gtr`)

## Overview
The `data_collection_gtr` pipeline is designed to fetch and preprocess data from the Gateway to Research (GtR) API. This pipeline is a component of a Kedro project, focusing on efficiently gathering and preparing GtR data for analysis and research within the project.

## Modules
### `nodes.py`
This module contains the key functions and classes for interacting with the GtR API. It includes:
- **Fetching GtR Data (`fetch_gtr_data`)**: Retrieves data from the GtR API based on specified parameters.
- **Preprocessing Data to DataFrame (`preprocess_data_to_df`)**: Converts the fetched raw data into a structured pandas DataFrame.
- **GtR Data Preprocessor Class**: A class that provides methods for preprocessing various types of GtR data, including organisations, funds, publications, and projects.

### `pipeline.py`
Defines the `data_collection_gtr` pipeline, which includes nodes for:
- Fetching data from various endpoints of the GtR API.
- Preprocessing the fetched data into a format suitable for further analysis and usage in the project.

### `utils.py`
Contains utility functions supporting the pipeline, including API configuration construction, main address extraction, and nested dictionary transformation.

## Usage
To execute the entire `data_collection_gtr` pipeline, run:
```bash
$ kedro run --pipeline data_collection_gtr
```

For running specific parts of the pipeline, such as a single endpoint, use the `--tags` flag. For example:
```bash
$ kedro run --pipeline data_collection_gtr --tags projects
```