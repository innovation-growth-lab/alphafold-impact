# Gateway to Research (GtR) Data Processing Pipeline (`data_processing_gtr`)

## Overview
The `data_processing_gtr` pipeline is part of a Kedro project aimed at processing data fetched from the Gateway to Research (GtR) API. This pipeline focuses on transforming raw GtR data into a structured format suitable for analysis and further use within the project.

## Modules
### `nodes.py`
This module contains the core functionalities for the pipeline:
- **Fetching GtR Data (`fetch_gtr_data`)**: Retrieves data from the GtR API using provided parameters.
- **Preprocessing Data to DataFrame (`preprocess_data_to_df`)**: Transforms the fetched data into a pandas DataFrame.
- **GtR Data Preprocessor Class**: A class providing methods for preprocessing various types of GtR data, including organisations, funds, publications, and projects.

### `pipeline.py`
Defines the `data_processing_gtr` pipeline, which includes:
- Nodes for fetching data from different GtR API endpoints.
- Nodes for preprocessing the fetched data, tailored to the specific structure of each data type from GtR.

## Usage
To execute the entire `data_processing_gtr` pipeline, run:
```bash
$ kedro run --pipeline data_processing_gtr
```

For running specific parts of the pipeline, such as a single data type, use the `--tags` flag. For example:
```bash
$ kedro run --pipeline data_processing_gtr --tags gtr.publications
```