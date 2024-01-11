from pathlib import Path
from typing import Dict, List

import pandas as pd

from src import df_col_to_ufloat, get_files_in_subdirectories, read_res_file


def process_files(res_dirs: Dict[str, List[Path]]) -> Dict[str, pd.DataFrame]:
    """
    method to process the res files from diffpy-CMI

    Args:
        res_dirs (Dict[str, List[str]]): dir containing res
        files from diffpy-CMI

    Returns:
        Dict[str, pd.DataFrame]: dict of dataframes with the results
    """
    results: Dict[str, pd.DataFrame] = {}
    for folder, files in res_dirs.items():
        dfs: List[pd.DataFrame] = []
        dfes: List[pd.DataFrame] = []
        if len(files) == 0:
            continue
        for file in files:
            df, dfe = read_res_file(Path(file))
            dfs.append(df)
            dfes.append(dfe)
        print(df)
        results[folder] = pd.concat(dfs)
        results[folder + "_errors"] = pd.concat(dfes)
        results[folder].reset_index(drop=True, inplace=True)
        results[folder + "_errors"].reset_index(drop=True, inplace=True)
    return results


def process_columns_and_save(
    res_dirs: Dict[str, List[Path]], results: Dict[str, pd.DataFrame]
):
    """
    method to combine the results and the errors.
    Using two dataframes and save them to a csv
    """
    for folder in res_dirs.keys():
        df = results[folder].copy()
        dfe = results[folder + "_errors"].copy()
        for col in df.columns:
            if col in ("name", "Rw"):
                continue
            #df = df_col_to_ufloat(df, dfe, col)
        df.to_csv(f"./table/{folder}.csv", index=False)
    return df


def main():
    """Creates tables from the res files including the errors"""
    directories = Path("./*/")
    extension = "res"
    res_dirs = get_files_in_subdirectories(directories, extension)
    print(res_dirs)
    results = process_files(res_dirs)
    process_columns_and_save(res_dirs, results)


if __name__ == "__main__":
    main()
