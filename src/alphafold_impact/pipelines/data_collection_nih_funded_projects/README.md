# NIH Funded Projects Data Collection Pipeline (`data_collection_nih_funded_projects`)

## Overview
The `data_collection_nih_funded_projects` Kedro pipeline collects NIH funded projects using the API at `https://api.reporter.nih.gov/v2/projects/search`.

## Modules
### `nodes.py`
This module contains functions for fetching and processing NIH funded projects data which are used in the pipeline. 

### `pipeline.py`
Defines the `data_collection_nih_funded_projects` pipeline, which fetches NIH funded projects data and saves it as CSVs in a partitioned dataset on S3.

## Usage
To execute the entire `data_collection_nih_funded_projects` pipeline, run:
```bash
$ kedro run --pipeline data_collection_nih_funded_projects
```