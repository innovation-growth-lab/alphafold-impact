# OpenAlex Data Collection Pipeline (`data_collection_oa`)

## Overview
The `data_collection_oa` pipeline is designed for fetching and preprocessing data from the OpenAlex (OA) API. This pipeline is part of a Kedro project and focuses on collecting data for research and analysis purposes.

## Modules
### `nodes.py`
This module contains the core functions for interacting with the OA API. Key functionalities include:
- **Collecting Papers (`collect_papers`)**: Fetches papers based on specified work IDs, with options for grouping, eager loading, and parallel processing.
- **Loading Work IDs (`load_work_ids`)**: Loads work IDs from a dataset and prepares them for further processing.
- **Retrieving Works for Concepts and Years (`retrieve_oa_works_for_concepts_and_years`)**: Gathers OA works based on concept IDs and publication years, handling API limitations.
- **Preprocessing DOI Values (`preprocess_publication_doi`)**: Prepares DOI values from Gateway to Research data for compatibility with OA filters.
- **Creating DOI Lists (`create_list_doi_inputs`)**: Generates lists of DOI values from preprocessed data.
- **Loading Referenced Work IDs (`load_referenced_work_ids`)**: Extracts work IDs and DOI mappings from a dataset.
- **Dynamically descending a citation tree (`fetch_citation_depth`)**: Iterates over a list of papers to yield OA responses.

### `pipeline.py`
Defines the `data_collection_oa` pipeline, which includes several sub-pipelines:
- **Base Pipeline**: Focuses on collecting OA publications for specific work IDs and their citations.
- **GtR Pipeline**: Targets OpenAlex publications matching DOIs from Gateway to Research and collects papers cited by these publications.
- **Concepts and Years Pipeline**: Gathers OA works for specific concepts and publication years.
- **Downstream Impact Pipeline**: Analyses the downstream impact of works based on citations and references.
- **Network pipeline**: Creates a network of paper citations from a seed paper.

### `utils.py`
Contains utility functions supporting data collection tasks, including parsing API responses, handling pagination, and error management.

## Usage
To execute the entire `data_collection_oa` pipeline, run:
```bash
$ kedro run --pipeline data_collection_oa
```

For running specific parts of the pipeline, use the `--namespace` and `--tags` flags. For example:
```bash
$ kedro run --pipeline data_collection_oa --namespace oa.data_collection.gtr
$ kedro run --pipeline data_collection_oa --tags downstream_impact
```
