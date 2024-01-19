# Data Collection GtR Pipeline

## Overview

This pipeline is designed to fetch and preprocess data from the Gateway to Research (GtR) API. It handles various data types such as organisations, funds, publications, and projects, efficiently processing them into a format suitable for in-depth analysis.

## Nodes Description

### fetch_gtr_data
- **Purpose**: Retrieves data (e.g., organisations, funds, publications, projects) from the GtR API based on specified parameters.
- **Inputs**: A dictionary of parameters and the endpoint URL.
- **Outputs**: A list of dictionaries containing the raw data fetched from the API.

### preprocess_data_to_df
- **Purpose**: Transforms the raw data into a structured pandas DataFrame and preprocesses it according to the data type.
- **Inputs**: Raw data list and the type of endpoint (e.g., organisations, funds).
- **Outputs**: A pandas DataFrame with preprocessed data.

### GtRDataPreprocessor Class
- **Purpose**: Provides methods for preprocessing different types of GtR data.
- **Methods**:
  - `_preprocess_organisations`: Processes organisation data, extracting main addresses and cleaning data.
  - `_preprocess_funds`: Processes funds data, extracting financial values and cleaning data.
  - `_preprocess_publications`: Processes publication data, extracting relevant information and restructuring.
  - `_preprocess_projects`: Processes project data, transforming nested dictionaries and cleaning data.

## Running the Pipeline

To execute the entire data collection pipeline, use:

```bash
$ kedro run --pipeline data_collection_gtr
```

For running specific parts of the pipeline, use the corresponding tags:

```bash
$ kedro run --tags organisations
$ kedro run --tags funds
$ kedro run --tags publications
$ kedro run --tags pronjects
```

This pipeline dynamically interacts with various GtR API endpoints as defined in the parameters, ensuring a comprehensive collection and processing of data.

---

**Note**: The `nodes.py` module contains the implementation details of the functions used in the pipeline. The `pipeline.py` script defines the structure and sequence of the data processing workflow. The `utils.py` module provides auxiliary functions that support data processing tasks.