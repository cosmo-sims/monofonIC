#!/usr/bin/env python3

import argparse
import configparser
import pathlib
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--monofonic", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("configs", nargs="+")
    arguments = parser.parse_args()

    workdir = pathlib.Path(arguments.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    generated_inputs = []
    processes = []

    for source_name in arguments.configs:
        source = pathlib.Path(source_name)
        config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
        config.optionxform = str
        config.read(source)
        config["cosmology"]["transfer"] = "CLASS"
        config["execution"]["WriteBackscalingDiagnostics"] = "false"
        destination = workdir / source.name
        with destination.open("w") as output:
            config.write(output)

        processes.append(subprocess.Popen(
            [arguments.monofonic, str(destination)],
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        ))
        generated_inputs.append(
            workdir / f"{destination.stem}_input_class_parameters.ini"
        )

    for process in processes:
        output, _ = process.communicate()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(
                process.returncode, process.args, output=output
            )

    reference = generated_inputs[0].read_bytes()
    for candidate in generated_inputs[1:]:
        if candidate.read_bytes() != reference:
            raise AssertionError(
                f"Target CLASS input differs between backscaling configurations: "
                f"{generated_inputs[0].name} and {candidate.name}"
            )

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (subprocess.CalledProcessError, AssertionError) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
