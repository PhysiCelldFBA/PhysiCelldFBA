"""
Quantify how degenerate the FBA flux distribution is in the cancer_core_metabolism
simulation, by replaying each cell's own LP offline and running FVA on it.

Why
---
Standard FBA has alternative optima: many flux distributions give the same growth
rate, and the solver returns whichever vertex the simplex lands on. In dFBA that
matters more than in a steady-state study, because the exchange fluxes are what
couple the cell to BioFVM. If the secretion profile is not unique, neither is the
extracellular environment it produces.

This module measures that freedom rather than assuming it. For every distinct set
of exchange bounds the simulation actually generated, it reconstructs the LP,
fixes the objective at its optimum, and asks how far each exchange flux could move
without changing growth. The width of that interval is the degeneracy.

How the LP is reconstructed
---------------------------
The addon constrains the problem in three places, and all three are reproduced
here or the replay would be answering a different question:

  1. each exchange reaction's lower bound is the Michaelis-Menten uptake limit
     computed from the local substrate concentration -- this is the per-cell,
     per-step quantity recorded into custom_data as "<flux_id>_lb";
  2. the objective reaction's upper bound is capped at max_growth_rate;
  3. the death model's trigger flux is given a lower bound equal to
     death_flux_threshold (here a maintenance ATP demand).

Exchange upper bounds are never touched by the simulation, so they keep their
static values from the SBML.

Identifiers keep their SBML prefixes (R_, M_), matching the addon, so cobrapy is
loaded with its identifier rewriting disabled.

Usage
-----
    python3 fva_degeneracy.py --output-dir ../output --sbml ../config/MCF7_Zielinski_2017_core_trimmed.sbml

Requires cobra, pandas, numpy, scipy. On a machine with several Python
installations, use the interpreter that can import cobra.
"""

from __future__ import annotations

import argparse
import glob
import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
import scipy.io as sio

# --------------------------------------------------------------------------- #
# Configuration of the simulation being replayed.
# These must match the <intracellular> block of the configuration that produced
# the output, or the reconstructed LP is not the one the cells actually solved.
# --------------------------------------------------------------------------- #

@dataclass
class ModelSetup:
    objective_reaction: str = "R_biomass_reaction"
    max_growth_rate: float = 0.026803          # objective reaction upper bound
    death_trigger_flux: Optional[str] = "R_ATPM"
    death_flux_threshold: float = 1.07         # lower bound on the trigger flux
    exchanges: Sequence[str] = field(default_factory=lambda: [
        "R_EX_arg_L_e", "R_EX_asn_L_e", "R_EX_asp_L_e", "R_EX_chol_e",
        "R_EX_co2_e", "R_EX_glc_e", "R_EX_gln_L_e", "R_EX_glu_L_e",
        "R_EX_gly_e", "R_EX_ile_L_e", "R_EX_lac_L_e", "R_EX_leu_L_e",
        "R_EX_lys_L_e", "R_EX_nh4_e", "R_EX_o2_e", "R_EX_orn_e",
        "R_EX_phe_L_e", "R_EX_pi_e", "R_EX_ser_L_e", "R_EX_thr_L_e",
        "R_EX_trp_L_e", "R_EX_tyr_L_e", "R_EX_val_L_e",
    ])


# --------------------------------------------------------------------------- #
# Reading PhysiCell cell output
# --------------------------------------------------------------------------- #

def expand_labels(labels_node: ET.Element) -> List[str]:
    """Expand MultiCellDS cell labels into one name per matrix row.

    A label with size > 1 covers that many consecutive rows. Those are indexed
    explicitly here as ``name[i]``.

    This deliberately does not use the ``_x``/``_y``/``_z`` suffix convention
    found in some readers: that list has only three entries, while several of
    these arrays are sized by the number of substrates (23 in this project).
    Past the third element such readers either raise IndexError or emit
    duplicate names, in which case only the last substrate's values survive and
    the label list no longer aligns with the matrix rows.
    """
    names: List[str] = []
    for label in labels_node.findall("label"):
        size = int(label.get("size"))
        name = (label.text or "").strip()
        if size == 1:
            names.append(name)
        else:
            names.extend(f"{name}[{i}]" for i in range(size))
    return names


def read_cells_snapshot(xml_path: str) -> pd.DataFrame:
    """Read one MultiCellDS snapshot into a DataFrame with one row per cell."""
    root = ET.parse(xml_path).getroot()
    node = root.find(".//cell_population/custom/simplified_data")
    if node is None:
        raise ValueError(f"{xml_path}: no simplified_data node")

    labels = expand_labels(node.find("labels"))
    mat_name = node.find("filename").text
    mat_path = os.path.join(os.path.dirname(xml_path), mat_name)
    data = sio.loadmat(mat_path)["cells"]

    if data.shape[0] != len(labels):
        raise ValueError(
            f"{xml_path}: {len(labels)} labels but {data.shape[0]} matrix rows. "
            "The label expansion and the matrix have diverged; do not trust any "
            "column past the first multi-element variable."
        )

    time_node = root.find(".//current_time")
    time = float(time_node.text) if time_node is not None else np.nan

    df = pd.DataFrame(data.T, columns=labels)
    df.insert(0, "time", time)

    # Readable aliases for the position triple, which expand_labels names
    # position[0..2]. x is depth into the tissue from the vessel in this project.
    for i, axis in enumerate("xyz"):
        col = f"position[{i}]"
        if col in df.columns:
            df[axis] = df[col]
    return df


def read_run(output_dir: str, every: int = 1) -> pd.DataFrame:
    """Read every snapshot in an output directory, concatenated.

    ``every`` subsamples snapshots (``every=6`` keeps one in six).
    """
    xml_files = sorted(glob.glob(os.path.join(output_dir, "output*.xml")))
    if not xml_files:
        raise FileNotFoundError(f"no output*.xml in {output_dir}")
    frames = [read_cells_snapshot(f) for f in xml_files[::every]]
    df = pd.concat(frames, ignore_index=True)
    print(f"read {len(frames)} snapshots from {output_dir}: "
          f"{len(df)} cell-timepoints, t = {df.time.min():.0f} to {df.time.max():.0f} min")
    return df


# --------------------------------------------------------------------------- #
# Deduplicating the recorded bound sets
# --------------------------------------------------------------------------- #

def bound_matrix(df: pd.DataFrame, exchanges: Sequence[str]) -> np.ndarray:
    missing = [f"{e}_lb" for e in exchanges if f"{e}_lb" not in df.columns]
    if missing:
        raise KeyError(
            f"missing recorded lower bounds: {missing}. The simulation must be run "
            "with the custom_data variables that record them."
        )
    return df[[f"{e}_lb" for e in exchanges]].to_numpy(dtype=float)


def unique_bound_sets(df: pd.DataFrame, exchanges: Sequence[str],
                      decimals: int = 10) -> tuple:
    """Collapse identical bound vectors.

    Returns ``(unique_bounds, index_per_row)`` where ``index_per_row[i]`` is the
    row of ``unique_bounds`` that cell-timepoint ``i`` maps to.

    The compression is usually large, because the Michaelis-Menten bound is a
    function of the local substrate concentration alone -- not of the cell -- so
    every cell sharing a voxel at a given time has an identical bound vector.
    The exception is an exchange whose uptake was clamped to the amount actually
    present in the voxel, which does depend on the cell's dry weight.

    ``decimals`` sets the rounding applied before comparison. The default is
    fine enough to be effectively exact; coarsen it to trade a little accuracy
    for fewer FVA problems, and report whatever value was used.
    """
    bounds = np.round(bound_matrix(df, exchanges), decimals=decimals)
    uniq, inverse = np.unique(bounds, axis=0, return_inverse=True)
    print(f"bound sets: {len(bounds)} recorded -> {len(uniq)} unique "
          f"(x{len(bounds) / max(len(uniq), 1):.1f} compression, rounded at 1e-{decimals})")
    return uniq, inverse


# --------------------------------------------------------------------------- #
# The replay
# --------------------------------------------------------------------------- #

def load_model(sbml_path: str):
    """Load the SBML with cobrapy, keeping the R_/M_ identifier prefixes."""
    from cobra.io import read_sbml_model
    noop = lambda x: x
    model = read_sbml_model(sbml_path, f_replace={
        "F_GENE": noop, "F_GENE_REV": noop,
        "F_REACTION": noop, "F_REACTION_REV": noop,
        "F_SPECIE": noop, "F_SPECIE_REV": noop,
    })
    print(f"loaded {sbml_path}: {len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites, {len(model.genes)} genes")
    return model


def apply_simulation_constraints(model, setup: ModelSetup) -> None:
    """Apply the constraints the addon imposes once, at initialisation.

    Call this on the model (not inside a per-bound-set context) so it applies to
    every replay.
    """
    obj = model.reactions.get_by_id(setup.objective_reaction)
    obj.upper_bound = setup.max_growth_rate
    model.objective = obj
    if setup.death_trigger_flux:
        trig = model.reactions.get_by_id(setup.death_trigger_flux)
        trig.lower_bound = setup.death_flux_threshold
        print(f"constraints: {obj.id}.ub = {obj.upper_bound}, "
              f"{trig.id}.lb = {trig.lower_bound}")


def fva_for_bound_sets(model, setup: ModelSetup, uniq_bounds: np.ndarray,
                       fraction_of_optimum: float = 1.0,
                       progress_every: int = 50) -> pd.DataFrame:
    """Run FVA over the exchange reactions for each unique bound set.

    For each set: apply the recorded lower bounds, maximise the objective to get
    Z*, then find the min and max of every exchange flux subject to the objective
    staying at (a fraction of) Z*. The width max - min is the degeneracy of that
    flux at that operating point.

    Infeasible bound sets are recorded with status "infeasible" and NaN ranges;
    they correspond to cells the simulation also found infeasible, i.e. starving.
    """
    from cobra.flux_analysis import flux_variability_analysis

    exchanges = list(setup.exchanges)
    records = []
    for k, lbs in enumerate(uniq_bounds):
        if progress_every and k % progress_every == 0:
            print(f"  FVA {k}/{len(uniq_bounds)}", flush=True)
        with model:
            for rid, lb in zip(exchanges, lbs):
                model.reactions.get_by_id(rid).lower_bound = float(lb)

            sol = model.optimize()
            if sol.status != "optimal" or sol.objective_value is None:
                records.append({"set": k, "status": sol.status, "growth": np.nan})
                continue
            growth = float(sol.objective_value)

            try:
                fva = flux_variability_analysis(
                    model, reaction_list=exchanges,
                    fraction_of_optimum=fraction_of_optimum, loopless=False)
            except Exception as exc:                      # solver trouble, not biology
                records.append({"set": k, "status": f"fva_failed: {exc}",
                                "growth": growth})
                continue

            row = {"set": k, "status": "optimal", "growth": growth}
            for rid in exchanges:
                lo = float(fva.at[rid, "minimum"])
                hi = float(fva.at[rid, "maximum"])
                row[f"{rid}_min"] = lo
                row[f"{rid}_max"] = hi
                row[f"{rid}_range"] = hi - lo
            records.append(row)

    out = pd.DataFrame.from_records(records).set_index("set")
    n_ok = int((out.status == "optimal").sum())
    print(f"FVA done: {n_ok}/{len(out)} bound sets optimal")
    return out


# --------------------------------------------------------------------------- #
# Self-check and summary
# --------------------------------------------------------------------------- #

def check_recorded_flux_within_range(df: pd.DataFrame, inverse: np.ndarray,
                                     fva: pd.DataFrame, setup: ModelSetup,
                                     tol: float = 1e-6) -> pd.DataFrame:
    """Verify that each flux the simulation recorded lies inside its FVA interval.

    This is the strongest available check that the replay reconstructed the same
    LP the simulation solved. A recorded flux outside its own feasible interval
    means the replay and the simulation disagree about the problem -- a missing
    constraint, a wrong bound, or an identifier mismatch -- not that FBA is wrong.
    """
    rows = []
    for rid in setup.exchanges:
        if f"{rid}_min" not in fva.columns:
            continue
        flux = df[rid].to_numpy(dtype=float)
        lo = fva[f"{rid}_min"].to_numpy()[inverse]
        hi = fva[f"{rid}_max"].to_numpy()[inverse]
        ok = np.isnan(lo) | ((flux >= lo - tol) & (flux <= hi + tol))
        rows.append({"reaction": rid,
                     "checked": int(np.sum(~np.isnan(lo))),
                     "outside": int(np.sum(~ok)),
                     "worst_excess": float(np.nanmax(
                         np.maximum(lo - flux, flux - hi))) if len(flux) else np.nan})
    out = pd.DataFrame(rows)
    bad = int(out.outside.sum())
    if bad:
        print(f"!! {bad} recorded fluxes fall outside their FVA interval -- the "
              "replay does not match the simulation's LP. Investigate before "
              "trusting any range below.")
    else:
        print("self-check passed: every recorded flux lies within its FVA interval")
    return out


def summarise(df: pd.DataFrame, inverse: np.ndarray, fva: pd.DataFrame,
              setup: ModelSetup) -> pd.DataFrame:
    """Per-exchange degeneracy summary over all cell-timepoints."""
    rows = []
    for rid in setup.exchanges:
        col = f"{rid}_range"
        if col not in fva.columns:
            continue
        rng = fva[col].to_numpy()[inverse]
        flux = np.abs(df[rid].to_numpy(dtype=float))
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(flux > 1e-12, rng / flux, np.nan)
        rows.append({
            "reaction": rid,
            "median_range": np.nanmedian(rng),
            "p95_range": np.nanpercentile(rng, 95) if np.any(~np.isnan(rng)) else np.nan,
            "max_range": np.nanmax(rng) if np.any(~np.isnan(rng)) else np.nan,
            "median_rel": np.nanmedian(rel),
            "frac_degenerate": float(np.nanmean(rng > 1e-9)),
        })
    out = pd.DataFrame(rows).sort_values("median_range", ascending=False)
    return out.reset_index(drop=True)


def attach_ranges(df: pd.DataFrame, inverse: np.ndarray, fva: pd.DataFrame,
                  setup: ModelSetup) -> pd.DataFrame:
    """Return df with the FVA range of each exchange joined onto every cell-timepoint."""
    out = df.copy()
    out["bound_set"] = inverse
    out["replay_status"] = fva["status"].to_numpy()[inverse]
    for rid in setup.exchanges:
        if f"{rid}_range" in fva.columns:
            out[f"{rid}_range"] = fva[f"{rid}_range"].to_numpy()[inverse]
            out[f"{rid}_fva_min"] = fva[f"{rid}_min"].to_numpy()[inverse]
            out[f"{rid}_fva_max"] = fva[f"{rid}_max"].to_numpy()[inverse]
    return out


# --------------------------------------------------------------------------- #

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output-dir", default="../output")
    ap.add_argument("--sbml", default="../config/MCF7_Zielinski_2017_core_trimmed.sbml")
    ap.add_argument("--every", type=int, default=1,
                    help="keep one snapshot in EVERY (default: all)")
    ap.add_argument("--decimals", type=int, default=10,
                    help="rounding applied before deduplicating bound sets")
    ap.add_argument("--fraction-of-optimum", type=float, default=1.0)
    ap.add_argument("--save-prefix", default="fva_degeneracy",
                    help="write <prefix>_cells.csv and <prefix>_summary.csv")
    args = ap.parse_args()

    setup = ModelSetup()
    df = read_run(args.output_dir, every=args.every)
    uniq, inverse = unique_bound_sets(df, setup.exchanges, decimals=args.decimals)

    model = load_model(args.sbml)
    apply_simulation_constraints(model, setup)
    fva = fva_for_bound_sets(model, setup, uniq,
                             fraction_of_optimum=args.fraction_of_optimum)

    check_recorded_flux_within_range(df, inverse, fva, setup)
    summary = summarise(df, inverse, fva, setup)
    print("\nper-exchange degeneracy (mmol/gDW/h), widest first:")
    print(summary.to_string(index=False))

    cells = attach_ranges(df, inverse, fva, setup)
    cells.to_csv(f"{args.save_prefix}_cells.csv", index=False)
    summary.to_csv(f"{args.save_prefix}_summary.csv", index=False)
    print(f"\nwrote {args.save_prefix}_cells.csv and {args.save_prefix}_summary.csv")


if __name__ == "__main__":
    main()
