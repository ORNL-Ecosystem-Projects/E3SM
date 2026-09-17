#!/usr/bin/env python3
"""Phase 2 parameter-file and closed-box utilities.

The closed-box implementation is an independent conservation oracle for the
ELM cascade contract: donor loss, pathway respiration, receiver gain, and the
corresponding mineral N/P residual. It is not a replacement for CIME tests.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any


POOLS = (
    "atmosphere",
    "litr1",
    "litr2",
    "litr3",
    "cwd",
    "soil1",
    "soil2",
    "soil3",
    "soil4",
    "bacteria",
    "fungi",
    "dom",
)

TRANSITIONS = (
    ("CWDL2", 4, 2), ("CWDL3", 4, 3),
    ("L1B", 1, 9), ("L1F", 1, 10), ("L1S1", 1, 5),
    ("L2B", 2, 9), ("L2F", 2, 10), ("L2S2", 2, 6),
    ("L3B", 3, 9), ("L3F", 3, 10), ("L3S3", 3, 7),
    ("S1B", 5, 9), ("S1F", 5, 10), ("S1S2", 5, 6),
    ("S2B", 6, 9), ("S2F", 6, 10), ("S2S3", 6, 7),
    ("S3B", 7, 9), ("S3F", 7, 10), ("S3S4", 7, 8),
    ("S4B", 8, 9), ("S4F", 8, 10),
    ("BS1", 9, 5), ("FS1", 10, 5),
    ("BS2", 9, 6), ("FS2", 10, 6),
    ("BS3", 9, 7), ("FS3", 10, 7),
    ("BS4", 9, 8), ("FS4", 10, 8),
    ("DOMS1", 11, 5), ("DOMS2", 11, 6),
    ("DOMS3", 11, 7), ("DOMS4", 11, 8),
    ("BDOM", 9, 11), ("FDOM", 10, 11),
    ("DOMB", 11, 9), ("DOMF", 11, 10),
    ("BATM", 9, 0), ("FATM", 10, 0),
    ("L1DOM", 1, 11), ("L2DOM", 2, 11), ("L3DOM", 3, 11),
    ("S1DOM", 5, 11), ("S2DOM", 6, 11),
    ("S3DOM", 7, 11), ("S4DOM", 8, 11),
)


def load_parameters(path: Path) -> tuple[dict[str, float], dict[str, Any]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    values = {
        name: float(metadata["value"])
        for scope in ("pft", "scalar")
        for name, metadata in document[scope].items()
    }
    return values, document


def parameter_provenance(document: dict[str, Any]) -> dict[str, str]:
    """Return a parameter-to-provenance-class mapping, rejecting ambiguity."""
    result: dict[str, str] = {}
    for provenance_class, metadata in document["provenance"].items():
        for name in metadata["parameters"]:
            if name in result:
                raise ValueError(f"duplicate provenance for {name}")
            result[name] = provenance_class
    return result


def cue_from_cn(substrate_cn: float, parameters: dict[str, float]) -> float:
    return parameters["CUEmax"] * min(
        1.0,
        parameters["microbe_cue_cn_target"]
        / (substrate_cn * parameters["CUEmax"]),
    )


def microbial_allocation(
    substrate_cn: float, parameters: dict[str, float]
) -> tuple[float, float]:
    exponent = parameters["microbe_allocation_cn_exponent"]
    bacteria = (parameters["cn_bacteria"] / substrate_cn) ** exponent
    fungi = (parameters["cn_fungi"] / substrate_cn) ** exponent
    bacteria_fraction = bacteria / (bacteria + fungi)
    return bacteria_fraction, 1.0 - bacteria_fraction


def build_pathways(
    parameters: dict[str, float], pool_cn: list[float]
) -> tuple[list[float], list[float]]:
    """Return path and respiration fractions in the Fortran transition order."""
    path = [0.0] * len(TRANSITIONS)
    respiration = [0.0] * len(TRANSITIONS)
    path[0:2] = [0.76, 0.24]

    def substrate_paths(
        indices: tuple[int, int, int, int],
        substrate_cn: float,
        dom_fraction: float,
        direct_fraction: float,
        cue: float,
    ) -> None:
        bacteria, fungi = microbial_allocation(substrate_cn, parameters)
        available = 1.0 - dom_fraction - direct_fraction
        ib, iff, direct, dom = indices
        path[ib] = available * bacteria
        path[iff] = available * fungi
        path[direct] = direct_fraction
        path[dom] = dom_fraction
        respiration[ib] = 1.0 - cue
        respiration[iff] = 1.0 - cue

    substrate_paths((2, 3, 4, 40), 90.0, parameters["l1dom_f"], parameters["l1s1_f"], cue_from_cn(90.0, parameters))
    substrate_paths((5, 6, 7, 41), 90.0, parameters["l2dom_f"], parameters["l2s2_f"], cue_from_cn(90.0, parameters))
    substrate_paths((8, 9, 10, 42), 90.0, parameters["l3dom_f"], parameters["l3s3_f"], cue_from_cn(90.0, parameters))
    substrate_paths((11, 12, 13, 43), pool_cn[5], parameters["s1dom_f"], parameters["s1s2_f"], parameters["m_rf_s1m"])
    substrate_paths((14, 15, 16, 44), pool_cn[6], parameters["s2dom_f"], parameters["s2s3_f"], parameters["m_rf_s2m"])
    substrate_paths((17, 18, 19, 45), pool_cn[7], parameters["s3dom_f"], parameters["s3s4_f"], parameters["m_rf_s3m"])

    bacteria, fungi = microbial_allocation(pool_cn[8], parameters)
    path[20] = bacteria * (1.0 - parameters["s4dom_f"])
    path[21] = fungi * (1.0 - parameters["s4dom_f"])
    path[46] = parameters["s4dom_f"]
    respiration[20] = respiration[21] = 1.0 - parameters["m_rf_s4m"]

    path[22] = parameters["m_bs1_f"]
    path[24] = parameters["m_bs2_f"]
    path[26] = parameters["m_bs3_f"]
    path[28] = 1.0 - sum(parameters[name] for name in ("m_batm_f", "m_bdom_f", "m_bs1_f", "m_bs2_f", "m_bs3_f"))
    path[34] = parameters["m_bdom_f"]
    path[38] = parameters["m_batm_f"]
    respiration[38] = 1.0

    path[23] = parameters["m_fs1_f"]
    path[25] = parameters["m_fs2_f"]
    path[27] = parameters["m_fs3_f"]
    path[29] = 1.0 - sum(parameters[name] for name in ("m_fatm_f", "m_fdom_f", "m_fs1_f", "m_fs2_f", "m_fs3_f"))
    path[35] = parameters["m_fdom_f"]
    path[39] = parameters["m_fatm_f"]
    respiration[39] = 1.0

    path[30] = parameters["m_doms1_f"]
    path[31] = parameters["m_doms2_f"]
    path[32] = parameters["m_doms3_f"]
    path[33] = 1.0 - sum(parameters[name] for name in ("m_domb_f", "m_domf_f", "m_doms1_f", "m_doms2_f", "m_doms3_f"))
    path[36] = parameters["m_domb_f"]
    path[37] = parameters["m_domf_f"]
    return path, respiration


def donor_sums(path: list[float]) -> dict[int, float]:
    sums = {pool: 0.0 for pool in range(1, len(POOLS))}
    for fraction, (_, donor, _) in zip(path, TRANSITIONS):
        sums[donor] += fraction
    return sums


def closed_box_step(
    carbon: list[float],
    nitrogen: list[float],
    phosphorus: list[float],
    mineral_n: float,
    mineral_p: float,
    respiration_total: float,
    rates: list[float],
    pool_cn: list[float],
    pool_cp: list[float],
    path: list[float],
    respiration: list[float],
) -> tuple[list[float], list[float], list[float], float, float, float]:
    dc = [0.0] * len(POOLS)
    dn = [0.0] * len(POOLS)
    dp = [0.0] * len(POOLS)

    for fraction, rf, (_, donor, receiver) in zip(path, respiration, TRANSITIONS):
        gross_c = carbon[donor] * rates[donor] * fraction
        gross_n = nitrogen[donor] * rates[donor] * fraction
        gross_p = phosphorus[donor] * rates[donor] * fraction
        received_c = gross_c * (1.0 - rf)

        dc[donor] -= gross_c
        dn[donor] -= gross_n
        dp[donor] -= gross_p
        respiration_total += gross_c * rf

        if receiver != 0:
            received_n = received_c / pool_cn[receiver]
            received_p = received_c / pool_cp[receiver]
            dc[receiver] += received_c
            dn[receiver] += received_n
            dp[receiver] += received_p
            mineral_n += gross_n - received_n
            mineral_p += gross_p - received_p
        else:
            mineral_n += gross_n
            mineral_p += gross_p

    carbon = [value + change for value, change in zip(carbon, dc)]
    nitrogen = [value + change for value, change in zip(nitrogen, dn)]
    phosphorus = [value + change for value, change in zip(phosphorus, dp)]
    return carbon, nitrogen, phosphorus, mineral_n, mineral_p, respiration_total


def inject_parameters(input_path: Path, output_path: Path, definition_path: Path) -> None:
    try:
        import netCDF4  # type: ignore
    except ImportError as error:
        raise SystemExit("netCDF4 is required; run this command in the ELM Docker image") from error

    _, document = load_parameters(definition_path)
    provenance = parameter_provenance(document)
    if input_path.resolve() != output_path.resolve():
        output_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(input_path, output_path)

    with netCDF4.Dataset(output_path, "a") as dataset:
        if "pft" not in dataset.dimensions:
            raise SystemExit("input parameter file has no pft dimension")
        scalar_dimensions = ("allpfts",) if "allpfts" in dataset.dimensions else ()
        for scope, dimensions in (("pft", ("pft",)), ("scalar", scalar_dimensions)):
            for name, metadata in document[scope].items():
                if name in dataset.variables:
                    variable = dataset.variables[name]
                    if variable.dimensions != dimensions:
                        raise SystemExit(f"{name} has dimensions {variable.dimensions}, expected {dimensions}")
                else:
                    variable = dataset.createVariable(name, "f8", dimensions)
                variable[...] = float(metadata["value"])
                variable.units = metadata["units"]
                variable.long_name = metadata["long_name"]
                variable.legacy_name = name
                variable.source_revision = document["source_revision"]
                variable.source_provenance = provenance[name]
                variable.science_status = "phase2_test_only"
        dataset.microbial_decomposition_schema = document["schema"]
        dataset.microbial_decomposition_warning = document["warning"]


def validate_parameter_file(path: Path, definition_path: Path) -> None:
    try:
        import netCDF4  # type: ignore
    except ImportError as error:
        raise SystemExit("netCDF4 is required; run this command in the ELM Docker image") from error

    _, document = load_parameters(definition_path)
    with netCDF4.Dataset(path) as dataset:
        missing = [
            name
            for scope in ("pft", "scalar")
            for name in document[scope]
            if name not in dataset.variables
        ]
        if missing:
            raise SystemExit("missing microbial parameters: " + ", ".join(missing))
        if getattr(dataset, "microbial_decomposition_schema", "") != document["schema"]:
            raise SystemExit("microbial decomposition schema attribute is missing or incompatible")


def main() -> None:
    default_definition = Path(__file__).with_name("phase2_test_parameters.json")
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    inject = subparsers.add_parser("inject-parameters")
    inject.add_argument("--input", type=Path, required=True)
    inject.add_argument("--output", type=Path, required=True)
    inject.add_argument("--definitions", type=Path, default=default_definition)

    validate = subparsers.add_parser("validate-parameter-file")
    validate.add_argument("--file", type=Path, required=True)
    validate.add_argument("--definitions", type=Path, default=default_definition)

    arguments = parser.parse_args()
    if arguments.command == "inject-parameters":
        inject_parameters(arguments.input, arguments.output, arguments.definitions)
    elif arguments.command == "validate-parameter-file":
        validate_parameter_file(arguments.file, arguments.definitions)


if __name__ == "__main__":
    main()
