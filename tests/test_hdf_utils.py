import hashlib
from pathlib import Path

import h5py
import numpy as np

ignored_header_entries = {"Git Branch", "Git Revision", "Git Tag"}


def header_dict_to_simple_dict(header_dict):
    meta_str = {}
    for key, value in header_dict.items():
        if key in ignored_header_entries:
            continue
        if isinstance(value, bytes):
            value = value.decode()
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if isinstance(value, np.float64) or isinstance(value, np.float32):
            value = float(value)
        if isinstance(value, np.int32) or isinstance(value, np.uint32) or isinstance(value, np.uint64):
            value = int(value)
        meta_str[key] = value
    return meta_str


def array_to_metadata(array: np.ndarray):
    arr_data = np.asarray(array)
    return {
        "shape": list(arr_data.shape),
        "dtype": arr_data.dtype.name,
        "hash": hashlib.sha256(arr_data.data).hexdigest(),
    }


def hdf5_to_metadata(input: Path):
    attrs = {}
    arrays = {}
    with h5py.File(input, 'r') as f:
        for col_name in f.keys():
            attrs[col_name] = header_dict_to_simple_dict(f[col_name].attrs)
            if isinstance(f[col_name], h5py.Group):
                for subcol_name in f[col_name].keys():
                    arr_name = col_name + "_" + subcol_name
                    arrays[arr_name] = array_to_metadata(f[col_name][subcol_name])
                    attrs[arr_name] = header_dict_to_simple_dict(f[col_name].attrs)
        else:
            arrays[col_name] = array_to_metadata(f[col_name])
    return {
        "attrs": attrs,
        "arrays": arrays
    }
