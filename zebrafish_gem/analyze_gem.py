"""Predict gene essentiality and substrate effects from a specified SBML GEM."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import cobra
from cobra.io import read_sbml_model
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

EPS = 1e-6
FRACTION = 0.01
CARBON_LIMIT = 60.0
INORGANIC_C = {"CO2", "CHO3", "CH2O3", "HCO3", "H2CO3", "CO3"}
ROOT = Path(__file__).resolve().parent


def solve(model):
    growth = model.slim_optimize()
    status = model.solver.status
    if status != "optimal" or not math.isfinite(growth):
        return status, None
    return status, float(growth)


def profile_exchanges(model):
    records = []
    baseline = model.medium
    for reaction in model.exchanges:
        if len(reaction.metabolites) != 1:
            raise ValueError(f"Exchange {reaction.id} must have one metabolite")
        metabolite, coefficient = next(iter(reaction.metabolites.items()))
        formula = metabolite.formula
        carbon_count = metabolite.elements.get("C", 0) if formula and metabolite.elements else None
        records.append({
            "exchange_id": reaction.id, "metabolite": metabolite.name,
            "metabolite_id": metabolite.id, "formula": formula,
            "carbon_atoms": carbon_count, "coefficient": coefficient,
            "organic_carbon": bool(carbon_count) and formula not in INORGANIC_C,
            "default_feed": baseline.get(reaction.id, 0.0),
        })
    return records


def gene_screen(model, baseline):
    records = []
    for gene in model.genes:
        with model:
            gene.knock_out()
            status, growth = solve(model)
        fraction = growth / baseline if growth is not None else None
        if fraction is None:
            call = "unresolved"
        elif fraction < FRACTION:
            call = "essential"
        else:
            call = "nonessential"
        records.append({"gene_id": gene.id, "gene_name": gene.name,
                        "status": status, "biomass": growth,
                        "fraction_of_baseline": fraction, "call": call})
    return pd.DataFrame(records)


def carbon_screen(model, exchanges, baseline):
    native = model.medium
    # Exclude unspecified formulas from the background, as they may supply carbon.
    background = {row["exchange_id"]: row["default_feed"] for row in exchanges
                  if row["carbon_atoms"] is not None and not row["organic_carbon"]
                  and row["default_feed"] > 0}
    sole, omission = [], []
    with model:
        model.medium = background
        control_status, control_growth = solve(model)
        if control_status != "optimal" or abs(control_growth) > EPS:
            raise ValueError(f"Organic-carbon-free control: {control_status}, {control_growth}")
        for row in exchanges:
            if not row["organic_carbon"]:
                continue
            reaction = model.reactions.get_by_id(row["exchange_id"])
            cap = CARBON_LIMIT / (row["carbon_atoms"] * abs(row["coefficient"]))
            with model:
                if row["coefficient"] < 0:
                    reaction.lower_bound = -cap
                else:
                    reaction.upper_bound = cap
                status, growth = solve(model)
            call = "unresolved" if growth is None else (
                "supports_biomass" if growth - control_growth > EPS else "no_biomass")
            sole.append({**row, "feed_cap": cap, "status": status,
                         "biomass": growth, "call": call})
    for row in exchanges:
        exchange = row["exchange_id"]
        if row["organic_carbon"] and exchange in native:
            with model:
                model.medium = {key: value for key, value in native.items() if key != exchange}
                status, growth = solve(model)
            omission.append({**row, "status": status, "biomass": growth,
                             "fraction_of_baseline": growth / baseline if growth is not None else None})
    return pd.DataFrame(sole), pd.DataFrame(omission), control_growth


def plot(genes, sole, omission, baseline, out):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), layout="constrained",
                             gridspec_kw={"width_ratios": [1, 1, 1.8]})
    for axis, data, values, labels, title in [
        (axes[0], genes, ["essential", "nonessential", "unresolved"],
         ["Essential", "Nonessential", "Unresolved"], "A  Single-gene deletions"),
        (axes[1], sole, ["supports_biomass", "no_biomass", "unresolved"],
         ["Supports\nbiomass", "No\nbiomass", "Unresolved"], "B  Sole organic-carbon screen"),
    ]:
        counts = [int((data.call == value).sum()) for value in values]
        bars = axis.bar(labels, counts, color=["#127C82", "#7C91A5", "#D19945"])
        axis.bar_label(bars, padding=3)
        axis.set_ylim(0, max(counts or [1]) * 1.2 or 1)
        axis.set_ylabel("Count")
        axis.set_title(title, loc="left", fontweight="bold")
        axis.tick_params(axis="x", labelsize=9)
    ax = axes[2]
    reduced = omission[omission.biomass < baseline - EPS]
    top = reduced.nsmallest(10, "fraction_of_baseline")
    labels = [f"{row.metabolite[:32]} ({row.exchange_id})" for row in top.itertuples()]
    values = 100 * top.fraction_of_baseline
    bars = ax.barh(labels, values, color="#127C82")
    ax.bar_label(bars, labels=[f"{value:.2f}%" for value in values], padding=4, fontsize=8)
    ax.axvline(100, lw=1, color="#7C91A5", linestyle="--")
    ax.set_xlim(0, 114)
    ax.invert_yaxis()
    ax.set_xlabel("Biomass retained (% of default medium)")
    ax.set_title("C  Largest single-substrate omission effects", loc="left", fontweight="bold")
    fig.suptitle("Zebrafish-GEM  |  Flux balance analysis", fontsize=15, fontweight="bold")
    fig.savefig(out / "prediction_summary.png", dpi=180)
    fig.savefig(out / "prediction_summary.svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=ROOT / "models/Zebrafish-GEM.xml")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    model = read_sbml_model(str(args.model))
    model.solver = "glpk"
    status, baseline = solve(model)
    if status != "optimal" or baseline <= EPS:
        raise ValueError(f"Default-medium baseline: {status}, {baseline}")
    with model:
        model.medium = {}
        closed_status, closed_growth = solve(model)
    if closed_status != "optimal" or abs(closed_growth) > EPS:
        raise ValueError(f"Uptake-closed control: {closed_status}, {closed_growth}")
    exchanges = profile_exchanges(model)
    print(f"{model.id}: baseline={baseline:.8g}, genes={len(model.genes)}", flush=True)
    genes = gene_screen(model, baseline)
    print("Gene screen complete", flush=True)
    sole, omission, carbon_control = carbon_screen(model, exchanges, baseline)
    args.out.mkdir(parents=True, exist_ok=True)
    for name, frame in {"gene_deletion": genes,
                        "essential_genes": genes[genes.call == "essential"],
                        "sole_organic_carbon": sole, "substrate_omission": omission,
                        "exchanges": pd.DataFrame(exchanges)}.items():
        frame.to_csv(args.out / f"{name}.csv", index=False)
    with args.model.open("rb") as model_file:
        sha256 = hashlib.file_digest(model_file, "sha256").hexdigest()
    report = {
        "model": model.id, "sha256": sha256,
        "cobra_version": cobra.__version__, "solver": model.solver.interface.__name__,
        "objective": str(model.objective.expression),
        "reactions": len(model.reactions), "metabolites": len(model.metabolites),
        "genes": len(model.genes), "exchanges": len(exchanges),
        "default_feed_count": len(model.medium), "default_biomass": baseline,
        "uptake_closed_biomass": closed_growth, "no_organic_carbon_biomass": carbon_control,
        "essential_fraction_threshold": FRACTION, "growth_tolerance": EPS,
        "carbon_uptake_limit": CARBON_LIMIT, "gene_calls": dict(Counter(genes.call)),
        "sole_organic_carbon_calls": dict(Counter(sole.call)),
        "omission_conditions": len(omission),
        "omission_statuses": dict(Counter(omission.status)),
        "omission_reduced_biomass": int((omission.biomass < baseline - EPS).sum()),
    }
    (args.out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    plot(genes, sole, omission, baseline, args.out)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
