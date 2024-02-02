"""
This module contains the functions used in the data processing pipeline for
the NIH clinical trials data.

Functions:
    trials_with_references: Loads NIH clinical trials from partitioned
        dataset, filter out trials that do not reference a paper,
        concatenate results into one DataFrame.
"""
from typing import Dict, Callable, Any
import logging
import pandas as pd

logger = logging.getLogger(__name__)


def trials_with_references(
    nih_clinical_trials_partitions: Dict[str, Callable[[], Any]]
) -> pd.DataFrame:
    """Loads NIH clinical trials from partitioned dataset, filter out trials
    that do not reference a paper, concatenate results into one DataFrame.

    Args:
        nih_clinical_trials_partitions (dict): Partitioned dataset of NIH clinical trials.

    Returns:
        DataFrame of NIH clinical trials with a paper reference.
    """
    ct_with_refs_combined = pd.DataFrame()
    n_partitions = len(nih_clinical_trials_partitions)

    for idx, (_, ct_partition_load_func) in enumerate(
        sorted(nih_clinical_trials_partitions.items()), start=1
    ):
        ct_partition_data_with_refs = ct_partition_load_func().query(
            "`ProtocolSection.ReferencesModule.ReferenceList.Reference`.notna()"
        )

        ct_with_refs_combined = pd.concat(
            [ct_with_refs_combined, ct_partition_data_with_refs],
            ignore_index=True,
            join="outer",
        )

        logger.info(
            "Concatenated clinical trials with references from partition %s/%s",
            idx,
            n_partitions,
        )

    return ct_with_refs_combined
