# Gateway to Research (GtR) Data Processing Pipeline (`data_processing_gtr`)

## Overview
The `data_processing_gtr` pipeline aims at processing data fetched from the Gateway to Research (GtR) API. This pipeline focuses on transforming raw GtR data into a structured format suitable for analysis and further use within the project.

## Modules
### `nodes.py`
This module contains the core functionalities for the pipeline:
- **Fetching GtR Data (`load_gtr_data`)**: Retrieves data from the GtR API using provided parameters.
- **Preprocessing Institutions data (`process_institutions`)**: Transforms the fetched institutional data.
- **Preprocessing Publications data (`process_publications`)**: Transforms the fetched publications data.
- **Preprocessing Publications data (`process_alphafold_citations`)**: Collects references to direct or indirect citations.

## Usage
To execute the entire `data_processing_gtr` pipeline, run:
```bash
$ kedro run --pipeline data_processing_gtr
```

For running specific parts of the pipeline, such as a single data type, use the `--tags` flag. For example:
```bash
$ kedro run --pipeline data_processing_gtr --tags gtr.publications
```