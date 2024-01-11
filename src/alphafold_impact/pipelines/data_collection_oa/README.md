# Data Collection OpenAlex Pipeline

## Overview
This pipeline collects data from the OpenAlex API.

This pipeline has multiple pipeline components to allow for:
   - collecting incoming (cites) or outgoing (references) citations for a given work using the OpenAlex API. The pipeline is structured to create two separate flows – `incoming_pipeline` and `outgoing_pipeline` – for collecting incoming and outgoing citations respectively. Additionally, it includes a `downstream_impact_pipeline` to assess the downstream impact of a work.
   - retrieving works from specified related concept IDs and publication years.

## Running the OpenAlex Incoming and Outgoing Citation Collection Pipeline
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

5. **Run for Retrieving Works from Concept IDs and Publication Years:**
   ```
   $ kedro run --pipeline data_collection_oa --tags=works_for_concepts_and_years
   ```

## Pipeline Structure
The main pipeline is a template that is used to generate the two specific pipelines. It utilizes the `modular_pipeline` function to create the `incoming_pipeline` and `outgoing_pipeline`. Each pipeline is tagged accordingly, and the output folders are named dynamically based on the direction of the citations.

### Nodes
The key functions in the `nodes` module include:

Functions:
- `collect_citation_papers`: Collects all papers cited by specific work IDs.
- `load_work_ids`: Loads the file corresponding to a particular work_id in a PartitionedDataSet,
        extracts all ids, and returns these as a list.
- `retrieve_oa_works_for_concepts_and_years`: Retrieves OpenAlex works for the specified
        concept IDs and publication years.

Internal functions:
- `_revert_abstract_index`: Reverts the abstract inverted index to the original text.
- `_parse_results`: Parses OpenAlex API response to retain basic variables.
- `_citation_works_generator`: Creates a generator that yields a list of works from the OpenAlex API based on a
        given work ID.
- `_create_concept_year_filter`: Creates an API query filter string for the OpenAlex API.
- `_retrieve_oa_works_chunk`: Retrieves a chunk of OpenAlex works for a chunk of concept IDs.
- `_chunk_list`: Divides the input list into chunks of specified size.

## Novelties
### Modular Pipelines
Modular pipelines are used to create flexible, reusable, and easily maintainable pipeline sections. This project utilizes modular pipelines to separate the logic for collecting incoming and outgoing citations.

### Partitioned Datasets
The pipeline leverages Kedro's Partitioned Datasets to manage data storage efficiently. This allows handling large datasets by dividing them into manageable parts, each corresponding to a specific work ID.

## To Do
- **Saving Downstream Impact Data:** A desired feature is to save downstream impact data in a subfolder named after the `work_id`. This is currently achievable using runtime parameters but not with static ones. Future iterations of the pipeline may explore custom configurations to enable this functionality.

- **Custom Configurations:** Investigate potential methods for creating a custom config that allows dynamic naming of output directories based on static parameters.

