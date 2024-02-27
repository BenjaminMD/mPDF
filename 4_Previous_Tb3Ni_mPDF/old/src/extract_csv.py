from pathlib import Path
from typing import Dict, Tuple, Union

import numpy as np
import pandas as pd


def read_res_file(file_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # skip first 14 rows and Read until "Fixed Variables"
    rw_raw = np.genfromtxt(
        file_path, delimiter="\t", skip_header=11, max_rows=1, dtype=str
    )
    lines = np.genfromtxt(
        file_path, delimiter="\t", skip_header=15, max_rows=40, dtype=str
    )
    idx = [i for i, _ in enumerate(lines) if _.startswith("Fixed")]
    lines = np.array(lines[: idx[0]])
    # mixed types dict
    param: Dict[str, Union[float, str]] = {}
    error: Dict[str, Union[float, str]] = {}
    name = file_path.stem
    param["name"] = name
    error["name"] = name

    rw_raw_str = str(rw_raw)
    param["Rw"] = float(rw_raw_str.split()[1])
    for line in lines:
        p, v, _, err = line.split()
        if p.find("Biso") != -1:
            ...
            # pUiso = p.replace("Biso", "Uiso")
            # param[pUiso] = float(v) / (8 * np.pi**2)
        param[p] = float(v)
        error[p] = float(err)
    dfp = pd.DataFrame(param, index=[0])
    dfe = pd.DataFrame(error, index=[0])
    return dfp, dfe


def combine_split_and_save(dfs, path, name):
    df = pd.concat(dfs, ignore_index=True)

    for config in ["d4_t_con", "d4_s_con", "d4config"]:
        df_select = df[df["name"].str.startswith(config)]
        df_select = df_select.dropna(axis=1, how="all")

        if df_select.empty:
            continue

        df_select.to_csv(f"{path}/{config}_{name}.csv", index=False)
