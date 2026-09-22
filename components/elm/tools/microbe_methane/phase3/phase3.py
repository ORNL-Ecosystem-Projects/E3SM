#!/usr/bin/env python3
"""Inject and validate Phase 3 revised-methane reference parameters."""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path
from typing import Any


def load_manifest(path: Path) -> dict[str, Any]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(document.get("parameters"), list):
        raise ValueError("phase3 manifest must contain a parameters list")
    names = [entry[0] for entry in document["parameters"]]
    if len(names) != len(set(names)):
        raise ValueError("phase3 manifest contains duplicate parameter names")
    return document


def parameter_map(document: dict[str, Any]) -> dict[str, tuple[Any, ...]]:
    return {entry[0]: tuple(entry[1:]) for entry in document["parameters"]}


def inject_parameters(input_path: Path, output_path: Path, manifest_path: Path) -> None:
    try:
        import netCDF4  # type: ignore
    except ImportError as error:
        raise SystemExit("netCDF4 is required; run this command in the ELM Docker image") from error

    if input_path.resolve() == output_path.resolve():
        raise SystemExit("input and output must differ; the standard parameter file is never edited in place")

    document = load_manifest(manifest_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_path, output_path)

    with netCDF4.Dataset(output_path, "a") as dataset:
        if "pft" not in dataset.dimensions:
            raise SystemExit("input parameter file has no pft dimension")
        dimensions = ("allpfts",) if "allpfts" in dataset.dimensions else ()
        for name, source_id, value, units, provenance in document["parameters"]:
            if name in dataset.variables:
                variable = dataset.variables[name]
                if variable.dimensions != dimensions:
                    raise SystemExit(
                        f"{name} has dimensions {variable.dimensions}, expected {dimensions}"
                    )
            else:
                variable = dataset.createVariable(name, "f8", dimensions)
            variable[...] = float(value)
            variable.units = units
            variable.long_name = name.removeprefix("microbe_methane_").replace("_", " ")
            variable.legacy_name = source_id
            variable.source_revision = document["source_revision"]
            variable.source_provenance = provenance
            variable.science_status = "phase3_reference_test_only"

        dataset.microbe_methane_schema = document["schema"]
        dataset.microbe_methane_warning = document["warning"]
        dataset.microbe_methane_source_revision = document["source_revision"]
        dataset.microbe_methane_source_parameter_file = document["source_parameter_file"]


def validate_parameter_file(path: Path, manifest_path: Path) -> None:
    try:
        import netCDF4  # type: ignore
    except ImportError as error:
        raise SystemExit("netCDF4 is required; run this command in the ELM Docker image") from error

    document = load_manifest(manifest_path)
    with netCDF4.Dataset(path) as dataset:
        missing = [entry[0] for entry in document["parameters"] if entry[0] not in dataset.variables]
        if missing:
            raise SystemExit("missing revised methane parameters: " + ", ".join(missing))
        if getattr(dataset, "microbe_methane_schema", "") != document["schema"]:
            raise SystemExit("revised methane schema attribute is missing or incompatible")
        expected_globals = {
            "microbe_methane_source_revision": document["source_revision"],
            "microbe_methane_source_parameter_file": document["source_parameter_file"],
            "microbe_methane_warning": document["warning"],
        }
        for attribute, expected_attribute in expected_globals.items():
            actual_attribute = getattr(dataset, attribute, None)
            if actual_attribute != expected_attribute:
                raise SystemExit(
                    f"global attribute {attribute} is {actual_attribute!r}, "
                    f"expected {expected_attribute!r}"
                )
        for name, source_id, expected, units, provenance in document["parameters"]:
            variable = dataset.variables[name]
            values = variable[...]
            if values.size != 1:
                raise SystemExit(f"{name} must contain exactly one scalar value")
            actual = float(values.flat[0])
            if not math.isclose(actual, float(expected), rel_tol=1.0e-13, abs_tol=0.0):
                raise SystemExit(f"{name} is {actual!r}, expected {expected!r}")
            expected_attributes = {
                "units": units,
                "legacy_name": source_id,
                "source_revision": document["source_revision"],
                "source_provenance": provenance,
                "science_status": "phase3_reference_test_only",
            }
            for attribute, expected_attribute in expected_attributes.items():
                actual_attribute = getattr(variable, attribute, None)
                if actual_attribute != expected_attribute:
                    raise SystemExit(
                        f"{name} attribute {attribute} is {actual_attribute!r}, "
                        f"expected {expected_attribute!r}"
                    )


def main() -> None:
    default_manifest = Path(__file__).with_name("phase3_reference_parameters.json")
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    inject = subparsers.add_parser("inject-parameters")
    inject.add_argument("--input", type=Path, required=True)
    inject.add_argument("--output", type=Path, required=True)
    inject.add_argument("--manifest", type=Path, default=default_manifest)

    validate = subparsers.add_parser("validate-parameter-file")
    validate.add_argument("--file", type=Path, required=True)
    validate.add_argument("--manifest", type=Path, default=default_manifest)

    arguments = parser.parse_args()
    if arguments.command == "inject-parameters":
        inject_parameters(arguments.input, arguments.output, arguments.manifest)
    else:
        validate_parameter_file(arguments.file, arguments.manifest)


if __name__ == "__main__":
    main()
