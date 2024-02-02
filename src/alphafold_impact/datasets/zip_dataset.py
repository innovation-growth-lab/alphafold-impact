"""A dataset that loads filenames from all files in a zip file."""
import io
import zipfile
from pathlib import PurePosixPath
from typing import Any, Dict, List
from kedro.io import AbstractDataset
from kedro.io.core import get_filepath_str, get_protocol_and_path
import fsspec


class ZipDataset(AbstractDataset):
    """A dataset that loads data from a zip file.

    Args:
        AbstractDataset: The abstract dataset class.
    """

    def __init__(
        self,
        filepath: str,
        credentials: Dict[str, Any] = None,
    ):
        """Creates a new instance of ``ZipDataset``.

        Args:
            filepath: The location of the image file to load / save data.
            credentials: Credentials required to get access to the underlying filesystem.
                E.g. for ``GCSFileSystem`` it should look like `{"token": None}`.
        """
        protocol, path = get_protocol_and_path(filepath)
        self._protocol = protocol
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol, **credentials)  # ,  **fs_args)

    def _load(self) -> List[str]:
        ids = []
        load_path = get_filepath_str(self._filepath, self._protocol)

        with self._fs.open(load_path, "rb") as f:
            with zipfile.ZipFile(io.BytesIO(f.read()), "r") as zip_ref:
                for filename in zip_ref.namelist():
                    if filename.endswith(".xml"):
                        ids.append(filename[:-4])
        return ids

    def _save(self, data: Any) -> None:
        raise NotImplementedError("Save method not implemented for ZipDataset")

    def _describe(self) -> Dict[str, Any]:
        return dict(filepath=self._filepath)
