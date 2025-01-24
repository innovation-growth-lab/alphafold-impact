import pyreadr
from kedro.io import AbstractDataset
import pandas as pd


class RdsDataset(AbstractDataset):
    """
    Custom Kedro dataset for reading RDS files using pyreadr.
    """

    def __init__(self, filepath: str):
        """
        Initializes the RdsDataset.

        Args:
            filepath: Path to the `.rds` file to be read.
        """
        self._filepath = filepath

    def _load(self) -> pd.DataFrame:
        """
        Loads data from the `.rds` file and returns it as a pandas DataFrame.

        Returns:
            A pandas DataFrame containing the data from the `.rds` file.
        """
        result = pyreadr.read_r(self._filepath)
        # Assuming the RDS file contains a single DataFrame object
        if len(result) == 1:
            return next(iter(result.values()))
        else:
            raise ValueError(
                "RDS file contains multiple objects. Unable to determine which to load."
            )

    def _save(self, data: pd.DataFrame) -> None:
        """
        Save functionality is not implemented for RdsDataset, as RDS files are read-only in this implementation.
        """
        raise NotImplementedError(
            "Saving to RDS files is not supported in this implementation."
        )

    def _describe(self) -> dict:
        """
        Describes the dataset.

        Returns:
            A dictionary containing the dataset description.
        """
        return {"filepath": self._filepath}
