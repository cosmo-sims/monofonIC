#!/usr/bin/env python3

import argparse
import configparser
import math
import pathlib
import subprocess
import sys

import numpy as np


def read_config(path):
    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    parser.optionxform = str.lower
    parser.read(path)
    return parser


def as_bool(value):
    return value.strip().lower() in {"1", "yes", "true", "on"}


def write_class_input(config, path, root):
    cosmology = config["cosmology"]
    h = float(cosmology["h0"]) / 100.0
    omega_m = float(cosmology["omega_m"])
    omega_b = float(cosmology["omega_b"])
    tcmb = float(cosmology.get("tcmb", "2.7255"))
    masses = [float(cosmology.get(f"m_nu{i}", "0")) for i in range(1, 4)]
    masses = [mass for mass in masses if mass > 1.0e-9]
    omega_nu = sum(masses) / (93.14 * h * h)
    use_radiation = as_bool(cosmology.get("backscalingradiation", "true"))
    use_neutrinos = as_bool(cosmology.get("backscalingneutrinos", "true"))
    clustering_cdm = omega_m - omega_b - (omega_nu if use_neutrinos else 0.0)
    n_ur = float(cosmology.get("n_ur", str(3.046 - len(masses))))

    lines = [
        f"root = {root}",
        "write background = yes",
        f"h = {h:.17g}",
        f"T_cmb = {tcmb:.17g}",
        f"Omega_b = {omega_b:.17g}",
        f"Omega_cdm = {clustering_cdm:.17g}",
        "Omega_k = 0",
        "Omega_Lambda = 0",
        f"w0_fld = {float(cosmology.get('w_0', '-1')):.17g}",
        f"wa_fld = {float(cosmology.get('w_a', '0')):.17g}",
        f"N_ur = {n_ur:.17g}",
    ]
    if use_neutrinos and masses:
        omega_species = [mass / (93.14 * h * h) for mass in masses]
        temperature = math.pow(4.0 / 11.0, 1.0 / 3.0)
        lines.extend([
            f"N_ncdm = {len(masses)}",
            "m_ncdm = " + ", ".join(f"{mass:.17g}" for mass in masses),
            "T_ncdm = " + ", ".join(f"{temperature:.17g}" for _ in masses),
            "Omega_ncdm = " + ", ".join(f"{omega:.17g}" for omega in omega_species),
        ])
    else:
        lines.append("N_ncdm = 0")
    path.write_text("\n".join(lines) + "\n")
    return len(masses) if use_neutrinos else 0


def interpolate_positive(a_class, values, a_sample):
    return np.exp(np.interp(np.log(a_sample), np.log(a_class), np.log(values)))


def compare_column(name, actual, expected, rtol, atol=1.0e-12):
    difference = np.abs(actual - expected)
    allowed = atol + rtol * np.abs(expected)
    failures = difference > allowed
    if np.any(failures):
        index = int(np.argmax(difference / allowed))
        raise AssertionError(
            f"{name}: {np.count_nonzero(failures)} samples exceed tolerance; "
            f"worst at row {index}: monofonIC={actual[index]:.16e}, "
            f"CLASS={expected[index]:.16e}, relerr="
            f"{difference[index] / max(abs(expected[index]), atol):.3e}"
        )


def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--monofonic", required=True)
    argument_parser.add_argument("--class", dest="class_executable", required=True)
    argument_parser.add_argument("--config", required=True)
    argument_parser.add_argument("--workdir", required=True)
    argument_parser.add_argument("--rtol", type=float, default=5.0e-5)
    arguments = argument_parser.parse_args()

    config_path = pathlib.Path(arguments.config).resolve()
    workdir = pathlib.Path(arguments.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    config = read_config(config_path)

    subprocess.run(
        [arguments.monofonic, str(config_path)],
        cwd=workdir,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    diagnostic = workdir / f"{config_path.stem}_backscaling.txt"
    monofonic = np.loadtxt(diagnostic)

    class_input = workdir / f"{config_path.stem}_class.ini"
    class_root = workdir / f"{config_path.stem}_class_"
    n_ncdm = write_class_input(config, class_input, class_root)
    subprocess.run(
        [arguments.class_executable, str(class_input)],
        cwd=workdir,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    background_files = list(workdir.glob(f"{class_root.name}*background.dat"))
    if not background_files:
        raise AssertionError("CLASS did not produce a background file")
    class_background = np.loadtxt(max(background_files, key=lambda path: path.stat().st_mtime))

    z_class = class_background[:, 0]
    a_class = 1.0 / (1.0 + z_class)
    h0 = class_background[-1, 3]
    rho_g = class_background[:, 8]
    rho_b = class_background[:, 9]
    rho_cdm = class_background[:, 10]
    rho_nu = np.zeros_like(rho_g)
    for species in range(n_ncdm):
        rho_nu += class_background[:, 11 + 2 * species]
    fld_index = 11 + 2 * n_ncdm
    rho_ur = class_background[:, fld_index + 2]

    a_sample = monofonic[:, 0]
    use_radiation = as_bool(config["cosmology"].get("backscalingradiation", "true"))
    matter = interpolate_positive(a_class, rho_b + rho_cdm, a_sample) / (h0 * h0)
    radiation = (
        interpolate_positive(a_class, rho_g + rho_ur, a_sample) / (h0 * h0)
        if use_radiation else np.zeros_like(a_sample)
    )
    neutrino = (
        interpolate_positive(a_class, rho_nu, a_sample) / (h0 * h0)
        if np.any(rho_nu > 0.0) else np.zeros_like(a_sample)
    )
    matter0 = (rho_b[-1] + rho_cdm[-1]) / (h0 * h0)
    radiation0 = (rho_g[-1] + rho_ur[-1]) / (h0 * h0) if use_radiation else 0.0
    neutrino0 = rho_nu[-1] / (h0 * h0)
    omega_de = 1.0 - matter0 - radiation0 - neutrino0
    w0 = float(config["cosmology"].get("w_0", "-1"))
    wa = float(config["cosmology"].get("w_a", "0"))
    dark_energy = (
        omega_de * np.power(a_sample, -3.0 * (1.0 + w0 + wa))
        * np.exp(-3.0 * wa * (1.0 - a_sample))
    )
    class_values = {
        "matter": matter,
        "radiation": radiation,
        "neutrino": neutrino,
        "dark_energy": dark_energy,
        "H/H0": np.sqrt(matter + radiation + neutrino + dark_energy),
    }
    columns = {"H/H0": 2, "matter": 3, "radiation": 4, "neutrino": 5, "dark_energy": 6}
    for name, column in columns.items():
        compare_column(name, monofonic[:, column], class_values[name], arguments.rtol)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (subprocess.CalledProcessError, AssertionError) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
