"""Custom JSONL Dataset for Kedro."""

import json
import logging
from pathlib import PurePosixPath
from typing import Any, Dict, List
from kedro.io import AbstractDataset
from kedro.io.core import get_filepath_str, get_protocol_and_path
import fsspec

logger = logging.getLogger(__name__)


class JSONLDataset(AbstractDataset):
    """A dataset class for loading JSONL (JSON Lines) files.

    Each line in the file should be a valid JSON object.
    """

    def __init__(
        self,
        filepath: str,
        credentials: Dict[str, Any] = None,
        load_args: Dict[str, Any] = None,
        save_args: Dict[str, Any] = None,
    ):
        """Creates a new instance of JSONLDataset.

        Args:
            filepath: The location of the JSONL file to load/save data.
            credentials: Credentials required to get access to the underlying filesystem.
            load_args: Additional arguments for loading.
            save_args: Additional arguments for saving.
        """
        protocol, path = get_protocol_and_path(filepath)
        self._protocol = protocol
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol, **(credentials or {}))
        self._load_args = load_args or {}
        self._save_args = save_args or {}

    def _load(self) -> List[Dict[str, Any]]:
        """Loads data from a JSONL file.

        Returns:
            List of dictionaries, one for each line in the JSONL file.
        """
        load_path = get_filepath_str(self._filepath, self._protocol)

        records = []
        with self._fs.open(load_path, mode="r", encoding="utf-8") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                try:
                    record = json.loads(line)
                    records.append(record)
                except json.JSONDecodeError as e:
                    logger.warning(
                        "Failed to parse JSON on line %d in %s: %s",
                        line_num,
                        load_path,
                        e,
                    )
                    continue

        logger.info("Loaded %d records from %s", len(records), load_path)
        return records

    def _save(self, data: List[Dict[str, Any]]) -> None:
        """Saves data to a JSONL file.

        Args:
            data: List of dictionaries to save.
        """
        save_path = get_filepath_str(self._filepath, self._protocol)

        with self._fs.open(save_path, mode="w", encoding="utf-8") as f:
            for record in data:
                json.dump(record, f, **self._save_args)
                f.write("\n")

        logger.info("Saved %d records to %s", len(data), save_path)

    def _describe(self) -> Dict[str, Any]:
        """Describes the dataset.

        Returns:
            Dictionary containing the dataset description.
        """
        return {
            "filepath": self._filepath,
            "protocol": self._protocol,
            "load_args": self._load_args,
            "save_args": self._save_args,
        }
