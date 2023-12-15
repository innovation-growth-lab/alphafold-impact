# Data Collection GtR Pipeline - README.md

## Overview

This pipeline is designed to fetch and preprocess data from the GtR API. It is structured into several nodes, each performing specific functions to transform the raw data into a structured format suitable for further analysis.

## Nodes

### fetch_gtr_data
- **Function**: Fetches data from the GtR API based on provided parameters.
- **Inputs**: Parameters dictionary.
- **Outputs**: List of dictionaries containing raw data.

### preprocess_data_to_df
- **Function**: Converts the raw data into a pandas DataFrame.
- **Inputs**: List of dictionaries (raw data).
- **Outputs**: Pandas DataFrame with raw data.

### preprocess_organisations
- **Function**: Processes organisations data by extracting main address and dropping unnecessary columns.
- **Inputs**: Pandas DataFrame of organisations data.
- **Outputs**: Processed pandas DataFrame.

## Running the Pipeline

To execute this pipeline, use the following command in your terminal:

```bash
$ kedro run --pipeline data_collection_gtr
```

The pipeline fetches data from the specified GtR API endpoints, preprocesses it into a pandas DataFrame, and applies specific transformations for organisation data, ensuring it's ready for downstream tasks.