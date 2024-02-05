# NIH Clinical Trials Data Processing Pipeline (`data_processing_nih_clinical_trials`)

## Overview
The `data_processing_nih_clinical_trials` pipeline processes raw data fetched from `clinicaltrials.gov`. This pipeline focuses on transforming the raw data into format suitable for further analysis within the project.

## Modules
### `nodes.py`
This module contains the core functionalities for the pipeline:
- `trials_with_references`: Loads NIH clinical trials from partitioned dataset, filters out trials that do not reference a paper, concatenates results into one DataFrame.
- `trials_links_to_papers`: Produces a DataFrame which can be used to link NIH clinical trials to research papers via PubMed ID or DOI. This function's expected input is the output from trials_with_references.


### `pipeline.py`
Defines the `data_processing_nih_clinical_trials` pipeline.

## Usage
To execute the entire `data_processing_nih_clinical_trials` pipeline, run:
```bash
$ kedro run --pipeline data_processing_nih_clinical_trials
```

For running specific parts of the pipeline, such as a single data type, use the `--tags` flag. For example:
```bash
$ kedro run --pipeline data_processing_nih_clinical_trials --tags trials_links_to_papers
```