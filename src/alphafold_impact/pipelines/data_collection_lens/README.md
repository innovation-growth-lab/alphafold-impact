# Lens API Data Collection Pipeline (`data_collection_lens`)

## Overview
The `data_collection_lens` pipeline fetches data from the Lens API. This pipeline outputs individual JSON files containing patents that match the permutation of years, months, and jurisdictions (EU and US). Due to limits with the Lens API, only 100k patents can be downloaded in a given month.

## Modules
### `nodes.py`
This module includes:

- **Create Request Form (`create_request_form`)**: Generates the JSON request body and headers for the Lens API, considering dynamic parameters such as month, year, and jurisdiction.
- **Fetch Lens Data (`fetch_lens_data`)**: Executes the data retrieval process from the Lens API utilising scroll requests. This function is responsible for managing API response pagination and handling rate limits.

### `pipeline.py`
The script assembles nodes in a logical sequence to:
- Create request forms for fetching data from various Lens API endpoints.
- Execute data fetching and handle API pagination, ensuring efficient data collection.

### Usage
To deploy the `data_collection_lens` pipeline, execute the following command:

```bash
$ kedro run --pipeline data_collection_lens
```