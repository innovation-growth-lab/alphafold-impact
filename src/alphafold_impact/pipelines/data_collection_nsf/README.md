**README.md for NSF Awards and Funding Opportunities Pipeline**

## Overview

This pipeline is tailored for extracting and processing data from the National Science Foundation (NSF) API. It focuses on gathering detailed information about NSF awards and funding opportunities, including archived funding opportunities. The goal is to provide a dataset that can be utilised for further analysis on AF's impact on NSF funded research and innovation.

## Pipeline Components

### Fetch NSF Award Data
- **Purpose**: Retrieves detailed information for specific NSF awards, including publication data associated with each award.
- **Outputs**: A dictionary containing detailed award information, ready for further analysis.

### Fetch NSF Data
- **Purpose**: Collects data for a given list of NSF awards for a specified award year, processing each award to extract relevant information.
- **Outputs**: A structured dataset of NSF awards, including key details such as award ID, funding amounts, and associated publications.

### Fetch NSF Archived Funding Opportunities
- **Purpose**: Gathers information on archived NSF funding opportunities from the NSF website, providing insights into historical funding trends and priorities.
- **Outputs**: A pandas DataFrame containing details of archived funding opportunities, including titles, synopses, and URLs.

## Running the Pipeline

To execute the entire NSF data collection pipeline, use:

```bash
$ kedro run --pipeline data_collection_nsf
```

For running specific parts of the pipeline, such as fetching data for a particular year or archived funding opportunities, use the corresponding tags:

```bash
$ kedro run --pipeline data_collection_nsf --tags 2018
$ kedro run --pipeline data_collection_nsf --tags archived_programmes
```