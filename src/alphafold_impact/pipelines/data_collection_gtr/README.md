# Data Collection GtR Pipeline - README.md

## Overview

This pipeline is designed for fetching and preprocessing data from the GtR API, handling multiple data types like organisations and funds. It efficiently processes the data into a format ready for further analysis.

## Nodes

### fetch_gtr_data
- **Function**: Fetches data (organisations, funds, etc.) based on parameters.
- **Inputs**: Parameter dictionary and endpoint specification.
- **Outputs**: List of dictionaries with raw data.

### preprocess_data_to_df
- **Function**: Converts raw data into a pandas DataFrame and preprocesses based on the data type.
- **Inputs**: Raw data list and endpoint type.
- **Outputs**: Preprocessed pandas DataFrame.

### preprocess_organisations (internal)
- **Function**: Processes organisations data, extracting main address and cleaning up.

### preprocess_funds (internal)
- **Function**: Processes funds data, extracting financial values and cleaning up.

## Running the Pipeline

Run the pipeline with:

```bash
$ kedro run --pipeline data_collection_gtr
```

You can run the smaller pipelines by using the pipeline tags:

```bash
$ kedro run --tags gtr_funds
```

This pipeline dynamically fetches and processes data from various GtR endpoints, as specified in the parameters.