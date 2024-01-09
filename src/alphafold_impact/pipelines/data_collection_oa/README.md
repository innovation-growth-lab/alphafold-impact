# OpenAlex Incoming and Outgoing Citation Collection Pipeline

## Overview
This pipeline is designed for collecting incoming (cites) or outgoing (references) citations for a given work using the OpenAlex API. The pipeline is structured to create two separate flows – `incoming_pipeline` and `outgoing_pipeline` – for collecting incoming and outgoing citations respectively. Additionally, it includes a `downstream_impact_pipeline` to assess the downstream impact of a work.

## Running the Pipeline
The pipeline can be executed in different modes based on the requirements:

1. **Standard Run:**
   ```
   $ kedro run --pipeline data_collection_oa
   ```

2. **Run with a Specific Work ID:**
   ```
   $ kedro run --pipeline data_collection_oa --params=get.work_id=<work_id>
   ```

3. **Run with a Specific Direction:**
   ```
   $ kedro run --pipeline data_collection_oa --tags=<direction>
   ```

4. **Run for Downstream Impact Analysis:**
   ```
   $ kedro run --pipeline data_collection_oa --tags=downstream_impact
   ```

## Pipeline Structure
The main pipeline is a template that is used to generate the two specific pipelines. It utilizes the `modular_pipeline` function to create the `incoming_pipeline` and `outgoing_pipeline`. Each pipeline is tagged accordingly, and the output folders are named dynamically based on the direction of the citations.

### Nodes
The key functions in the `nodes` module include:

- `_revert_abstract_index`: Reverts the abstract inverted index to original text.
- `_parse_results`: Parses OpenAlex API response to extract specific information.
- `_citation_works_generator`: Yields a list of works based on a given work ID.
- `collect_citation_papers`: Collects all papers cited by specific work IDs.
- `load_work_ids`: Loads and extracts work IDs from a specified dataset.

## Novelties
### Modular Pipelines
Modular pipelines are used to create flexible, reusable, and easily maintainable pipeline sections. This project utilizes modular pipelines to separate the logic for collecting incoming and outgoing citations.

### Partitioned Datasets
The pipeline leverages Kedro's Partitioned Datasets to manage data storage efficiently. This allows handling large datasets by dividing them into manageable parts, each corresponding to a specific work ID.

## To Do
- **Saving Downstream Impact Data:** A desired feature is to save downstream impact data in a subfolder named after the `work_id`. This is currently achievable using runtime parameters but not with static ones. Future iterations of the pipeline may explore custom configurations to enable this functionality.

- **Custom Configurations:** Investigate potential methods for creating a custom config that allows dynamic naming of output directories based on static parameters.