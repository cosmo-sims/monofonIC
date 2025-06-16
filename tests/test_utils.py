import configparser
import json
import shutil
import subprocess
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from pytest_check import check

from test_hdf_utils import hdf5_to_metadata

this_dir = Path(__file__).parent
monofonic_bin = this_dir.parent / "build" / "monofonIC"
configs_dir = this_dir / "configs"


def read_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    return config


class MonofonicRun:
    def __init__(self, config_name: str, output_file: str = None):
        self.config_name = config_name
        self.config_basename = config_name.split(".")[0] + "_"
        if output_file is None:
            self.config = read_config(configs_dir / self.config_name)
            self.output_file = configs_dir / self.config["output"]["filename"]
        else:
            self.output_file = configs_dir / output_file

    @property
    def log_file(self):
        return configs_dir / (self.config_basename + "log.txt")

    @property
    def input_powerspec_file(self):
        return configs_dir / (self.config_basename + "input_powerspec.txt")

    def input_powerspec(self):
        return np.loadtxt(self.input_powerspec_file)

    @property
    def input_transfer_file(self):
        return configs_dir / (self.config_basename + "input_transfer.txt")

    def input_transfer(self):
        return np.loadtxt(self.input_transfer_file)

    def output_file_metadata(self):
        return hdf5_to_metadata(self.output_file)

    def output_file_metadata_from_json(self):
        with (configs_dir / (self.config_basename + "output_meta.json")).open() as f:
            return json.load(f)

    @property
    def input_class_parameters_file(self):
        return configs_dir / (self.config_basename + "input_class_parameters.ini")

    def input_class_parameters(self):
        params = {}
        with self.input_class_parameters_file.open() as f:
            for line in f:
                key, value = line.split("=")
                key = key.strip()
                value = value.strip()
                params[key] = value
        return params

    def save_as_reference(self):
        for file in [self.log_file, self.input_powerspec_file, self.input_transfer_file,
                     self.input_class_parameters_file]:
            dir = file.parent
            name = file.name.split(self.config_basename)[-1]
            ref_file = dir / (self.config_basename + "ref_" + name)
            print(file.name, "->", ref_file.name)
            shutil.copyfile(file, ref_file)
        meta_outfile = configs_dir / (self.config_basename + "ref_output_meta.json")
        print(f"saving output info to {meta_outfile.name}")
        with meta_outfile.open("w") as f:
            json.dump(self.output_file_metadata(), f, indent=2, ensure_ascii=False)

    def run(self):
        out = subprocess.run([
            monofonic_bin,
            self.config_name,
        ], cwd=this_dir / "configs", capture_output=True)
        assert out.returncode == 0
        return out

    def get_reference(self):
        return MonofonicRun(self.config_basename + "ref.conf", self.output_file.name)


def compare_per_column(arr1: np.ndarray, arr2: np.ndarray, col_names: list[str]):
    for i, name in enumerate(col_names):
        np.testing.assert_array_equal(arr1[i], arr2[i], err_msg=name)


def compare_runs(run: MonofonicRun):
    reference = run.get_reference()
    with check:
        assert reference.input_class_parameters() == run.input_class_parameters()
    with check:
        compare_per_column(
            reference.input_transfer(), run.input_transfer(),
            ["k", "delta_c", "delta_b", "delta_m", "delta_bc"]
        )
    with check:
        compare_per_column(
            reference.input_powerspec(), run.input_powerspec(),
            ["k", "P_dtot", "P_dcdm", "P_dbar"]
        )
    ref_output_meta = reference.output_file_metadata_from_json()
    output_meta = run.output_file_metadata()
    with check:
        assert output_meta["attrs"] == ref_output_meta["attrs"]
    with check:
        assert output_meta["arrays"] == ref_output_meta["arrays"]


def plot_powerspectra(run: MonofonicRun):
    reference = run.get_reference()
    ref_power = reference.input_powerspec()
    run_power = run.input_powerspec()
    fig, axs = plt.subplots(2, sharex=True)
    ax_abs: Axes = axs[0]
    ax_rel: Axes = axs[1]
    ax_abs.plot(ref_power[::, 0], ref_power[::, 1], label="Reference")
    ax_abs.plot(run_power[::, 0], run_power[::, 1], label="Run")
    if np.allclose(run_power[::, 0], ref_power[::, 0]):
        ax_rel.plot(ref_power[::, 0], run_power[::, 1] / ref_power[::, 1])
    for ax in axs:
        ax.set_xscale("log")
    ax_abs.set_yscale("log")
    plt.show()
