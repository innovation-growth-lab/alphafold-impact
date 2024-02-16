**README.md for the AlphaFold Citation Analysis Pipeline**

## Overview

This pipeline is designed to analyse how the scientific community cites the AlphaFold paper, using Semantic Scholar to label citations as being background, method, or result-oriented. Due to limits in the API, the citation details for AlphaFold are extracted from the reference list of every paper that cites it, whereas for any subsequent papers we use the citations' details endpoint (this returns up to 10,000 results).

## Pipeline Components

### Load AlphaFold Citation IDs
- **Purpose**: Retrieves the list of IDs for papers citing AlphaFold, using data from open-access sources.
- **Outputs**: A list of citation IDs, including identifiers like DOI and PMID, for further analysis.

### Fetch Citation Details
- **Purpose**: Gathers detailed information on citations for specific papers, focusing on how AlphaFold's findings are being applied in later research.
- **Outputs**: Detailed records of each citation, including the citing paper's focus and how it relates to AlphaFold.

### Get AlphaFold Citation Details
- **Purpose**: Filters through the citation details to focus on those specifically citing AlphaFold, analyzing the context and intent behind these citations.
- **Outputs**: A structured overview of AlphaFold's direct and indirect influence on subsequent research, presented in a clear and accessible format.

## Running the Pipeline

To start analyzing AlphaFold's citation network, run:

```bash
$ kedro run --pipeline data_collection_s2
```

For more focused analyses, such as examining only direct citations or exploring specific research areas influenced by AlphaFold, use tags to run subsets of the pipeline:

```bash
$ kedro run --tags direct_citations
$ kedro run --tags research_impact
```