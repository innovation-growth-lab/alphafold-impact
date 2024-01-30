# NIH Clinical Trials Data Collection Pipeline (`data_collection_nih_clinical_trials`)

## Overview
The `data_collection_nih_clinical_trials` Kedro pipeline is designed to do a full download of the latest available data from `clinicaltrials.gov`.

## Modules
### `nodes.py`
This module contains functions for fetching and processing NIH clinical trials data which are used in the pipeline. 

### `pipeline.py`
Defines the `data_collection_nih_clinical_trials` pipeline, which fetches NIH clinical trials data and saves it as CSVs in a partitioned dataset on S3.

## Usage
To execute the entire `data_collection_nih_clinical_trials` pipeline, run:
```bash
$ kedro run --pipeline data_collection_nih_clinical_trials
```