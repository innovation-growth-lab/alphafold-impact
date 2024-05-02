"""
This is a boilerplate pipeline 'data_collection_pdb'
generated using Kedro 0.19.1
"""

import logging
from typing import Dict, List, Union
import requests
from requests.adapters import HTTPAdapter, Retry
import pandas as pd

logger = logging.getLogger(__name__)


def _fetch_pbd_ids() -> List[str]:
    """Fetch the list of PDB ids to be used for data collection"""
    # Fetch the list of PDB ids
    response = requests.get(
        "https://data.rcsb.org/rest/v1/holdings/current/entry_ids", timeout=10
    )
    return response.json()


def _query(ids: list) -> str:
    ids_string = ", ".join(f'"{id}"' for id in ids)
    return f"""
    {{
        entries(entry_ids: [{ids_string}]) {{
            rcsb_id
            rcsb_accession_info {{
            initial_release_date
            }}
            audit_author {{
            name
            }}
            rcsb_primary_citation {{
            pdbx_database_id_PubMed
            pdbx_database_id_DOI
            }}
        }}
    }}
    """


def fetch_pbd_details(config: Dict[str, str]) -> Dict[str, Union[str, List[str]]]:
    """Fetch the details of a PDB id"""

    pbd_ids = _fetch_pbd_ids()
    session = requests.Session()
    retries = Retry(
        total=config["max_retries"], backoff_factor=config["backoff_factor"]
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))
    url = "https://data.rcsb.org/graphql"

    # create sublists of 50 items
    pbd_ids = [pbd_ids[i : i + 50] for i in range(0, len(pbd_ids), 50)]

    outputs = []
    for ids in pbd_ids:
        logger.info("Fetching details for %s", ids)
        query = _query(ids)
        response = requests.get(url, params={"query": query}, timeout=10)
        data = response.json()
        data = pd.DataFrame(data["data"]["entries"])

        # get date
        data["rcsb_accession_info"] = data["rcsb_accession_info"].apply(
            lambda x: x["initial_release_date"]
        )

        # get an author list from the list of dictionaries with key "name"
        data["audit_authors"] = data["audit_author"].apply(
            lambda x: ", ".join([d.get("name", "") for d in x])
        )

        # get pmid
        data["pmid"] = data["rcsb_primary_citation"].apply(
            lambda x: x.get("pdbx_database_id_PubMed", "") if x is not None else ""
        )

        # get doi
        data["doi"] = data["rcsb_primary_citation"].apply(
            lambda x: x.get("pdbx_database_id_DOI", "") if x is not None else ""
        )

        # append
        outputs.append(data)

    # concatenate
    outputs = pd.concat(outputs)

    # make sure ids are strings
    outputs["rcsb_id"] = outputs["rcsb_id"].astype(str)
    outputs["pmid"] = outputs["pmid"].astype(str)
    outputs["doi"] = outputs["doi"].astype(str)

    return outputs
