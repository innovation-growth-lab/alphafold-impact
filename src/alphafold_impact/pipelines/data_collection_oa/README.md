# Data Collection OpenAlex Pipeline - README.md

## Overview

This pipeline collects data from the OpenAlex API.

## Running the Pipeline

Run the pipeline with:

```bash
$ kedro run --pipeline data_collection_oa
```

The pipeline creates the following data catalog items:
    * `works_for_concepts`