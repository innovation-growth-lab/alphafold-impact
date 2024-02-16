# NIH Funding Opportunities Data Collection Pipeline (`data_collection_nih_funding_opportunities`)

## Overview
The `data_collection_nih_funding_opportunities` Kedro pipeline collects NIH funding opportunities using the API at `https://search.grants.nih.gov/guide/api/`.

## Modules
### `nodes.py`
This module contains functions for fetching and processing NIH funding opportunities data which are used in the pipeline. 

### `pipeline.py`
Defines the `data_collection_nih_funding_opportunities` pipeline, which fetches NIH funding opportunities data and saves it as a CSV on S3.

## Usage
To execute the entire `data_collection_nih_funding_opportunities` pipeline, run:
```bash
$ kedro run --pipeline data_collection_nih_funding_opportunities
```