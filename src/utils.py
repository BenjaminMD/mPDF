""""Utility functions for the project."""
from glob import glob
from pathlib import Path
from typing import Dict, List

from uncertainties import ufloat  # type: ignore # pylint: disable=import-error


def get_files_in_subdirectories(
    directory: Path, file_extension: str
) -> Dict[str, List[Path]]:
    """
    Get dict of files in subdirectories of a directory,
    with a given extension file_extension.
    """

    dir_files: Dict[str, List[Path]] = {}
    for folder in glob(f"{directory}/"):
        folder_path: Path = Path(folder)
        folder_name: str = folder_path.name

        search_str: str = f"{folder_path}/**/*.{file_extension}"
        file_list: List[str] = glob(search_str, recursive=True)
        file_paths: List[Path] = [Path(file) for file in file_list]

        dir_files[folder_name] = file_paths
    return dir_files


def df_col_to_ufloat(df, dfe, col):
    """Convert a column of a dataframe to ufloats based on the error column."""
    uval = [ufloat(val, err) for val, err in zip(df[col], dfe[col])]
    uval_str = [f"{val:.1ufS}" for val in uval]
    df[col] = uval_str
    return df
