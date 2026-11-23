"""Compare published B. subtilis 168 GEMs under matched aerobic FBA conditions."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import cobra
from cobra.io import load_json_model, read_sbml_model, write_sbml_model
from cobra.util.solver import linear_reaction_coefficients
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
MODELS = ("gapseq", "AlphaGEM")
EPS = 1e-6
ESSENTIAL_FRACTION = 0.01
CARBON_LIMIT = 60.0
COLORS = {"gapseq": "#326289", "AlphaGEM": "#158476"}


def solve(model):
    value = model.slim_optimize()
    status = model.solver.status
    return status, float(value) if status == "optimal" and math.isfinite(value) else None


def exchange_index(model):
    return {r.id: r for r in model.exchanges}


def read_medium(table, label, exchanges):
    medium = {}
    for row in table.to_dict("records"):
        exchange_id = row[label]
        if not exchange_id:
            continue
        reaction = exchanges[exchange_id]
        compound, coefficient = next(iter(reaction.metabolites.items()))
        if len(reaction.metabolites) != 1 or not compound.elements:
            raise ValueError(f"Cannot validate medium exchange {exchange_id}")
        if compound.elements.get("C", 0) and compound.formula != "CO2":
            raise ValueError(f"Organic carbon in background: {exchange_id}")
        medium[exchange_id] = float(row["max_uptake"]) / abs(coefficient)
    return medium


def inputs_for_substrate(row, label, exchanges, medium):
    reaction = exchanges[row[label]]
    compound, coefficient = next(iter(reaction.metabolites.items()))
    if len(reaction.metabolites) != 1 or compound.elements.get("C") != row["carbon_atoms"]:
        raise ValueError(f"Carbon count mismatch for {reaction.id}")
    cap = CARBON_LIMIT / float(row["carbon_atoms"])
    return medium | {reaction.id: cap / abs(coefficient)}, cap


def assay_controls(model, medium):
    records = []
    for label, limits in (
        ("native", model.medium), ("closed", {}), ("no_organic_carbon", medium)
    ):
        with model:
            model.medium = limits
            status, growth = solve(model)
        if label != "native" and not (
            status == "infeasible" or (status == "optimal" and abs(growth) <= EPS)
        ):
            raise ValueError(f"{model.id} failed control {label}: {status}, {growth}")
        records.append({"condition": label, "status": status, "growth": growth})
    return records


def substrate_screen(model, label, menu, exchanges, medium):
    records = []
    for row in menu.to_dict("records"):
        conditions, cap = inputs_for_substrate(row, label, exchanges, medium)
        with model:
            model.medium = conditions
            status, growth = solve(model)
        if status == "infeasible":
            call = "infeasible"
        elif growth is None:
            call = "unresolved"
        else:
            call = "supports_growth" if growth > EPS else "zero_growth"
        reaction = exchanges[row[label]]
        compound = next(iter(reaction.metabolites))
        neighbors = sorted(r.id for r in compound.reactions if not r.boundary)
        records.append({"model": label, "substrate": row["substrate"],
                        "exchange": row[label], "extracellular_metabolite": compound.id,
                        "nonboundary_neighbors": "|".join(neighbors),
                        "carbon_atoms": row["carbon_atoms"],
                        "max_uptake": cap, "status": status, "growth": growth, "call": call})
    return pd.DataFrame(records)


def gene_accession(gene, label):
    if label == "AlphaGEM":
        return gene.id
    match = re.fullmatch(r"(?:sp|tr)_([A-Z0-9]+)_\1", gene.id)
    return match[1] if match else ""


def gene_screen(model, label, baseline, proteome):
    records = []
    for gene in model.genes:
        accession = gene_accession(gene, label)
        with model:
            gene.knock_out()
            status, growth = solve(model)
        fraction = growth / baseline if growth is not None else None
        if fraction is not None:
            call = "essential" if fraction < ESSENTIAL_FRACTION else "nonessential"
        elif status == "infeasible":
            call = "essential"
        else:
            call = "unresolved"
        records.append({"model": label, "model_gene_id": gene.id,
                        "uniprot": accession, "in_reference_proteome": accession in proteome,
                        "gene_names": proteome.get(accession, ""), "status": status,
                        "growth": growth, "fraction_of_baseline": fraction, "call": call})
    return pd.DataFrame(records)


def gene_comparison(frames):
    columns = ["uniprot", "model_gene_id", "gene_names", "call", "status"]
    indexes = []
    for label in MODELS:
        frame = frames[label]
        valid = frame[frame.in_reference_proteome & ~frame.uniprot.duplicated(keep=False)]
        renamed = valid[columns].set_index("uniprot").add_prefix(f"{label}_")
        indexes.append(renamed)
    return indexes[0].join(indexes[1], how="outer").sort_index()


def reaction_annotations(model):
    """Collect explicit BiGG labels on internal reactions, preserving one-to-many links."""
    boundary = {r.id for r in model.boundary}
    annotations = {}
    for reaction in model.reactions:
        if reaction.id in boundary:
            continue
        values = reaction.annotation.get("bigg.reaction", [])
        if isinstance(values, str):
            values = [values]
        for value in values:
            annotations.setdefault(value, set()).add(reaction.id)
    return annotations


def topology_summary(model, label, path, baseline):
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {
        "model": label, "original_file": path.name, "sha256": digest,
        "reactions": len(model.reactions), "metabolites": len(model.metabolites),
        "genes": len(model.genes), "exchange_reactions": len(model.exchanges),
        "gpr_reactions": sum(bool(r.gene_reaction_rule) for r in model.reactions),
        "objective": ";".join(r.id for r in linear_reaction_coefficients(model)),
        "atpm_lower_bound": model.reactions.ATPM.lower_bound if "ATPM" in model.reactions else None,
        "glucose_baseline": baseline,
    }


def plot(genes, shared, substrates, out):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 7.6),
                                     gridspec_kw={"width_ratios": [1, 2]}, layout="constrained")
    categories = ["Both essential", "gapseq only", "AlphaGEM only", "Neither essential"]
    a, b = (shared[f"{label}_call"] for label in MODELS)
    common = a.isin(["essential", "nonessential"]) & b.isin(["essential", "nonessential"])
    a, b = a[common].eq("essential"), b[common].eq("essential")
    counts = [int(x.sum()) for x in (a & b, a & ~b, ~a & b, ~a & ~b)]
    bars = left.barh(categories, counts, color=["#474B65", *COLORS.values(), "#AEB8BF"])
    left.bar_label(bars, padding=5)
    left.invert_yaxis()
    left.set_ylim(len(categories) + 0.4, -0.65)
    left.set_xlim(0, max(counts) * 1.20)
    left.set_xlabel("Count among shared UniProt genes")
    left.set_title(f"A  Essentiality agreement (n = {len(a):,})", loc="left", fontweight="bold")
    counts_text = "\n".join(
        f"{label}  {int((genes[label].call == 'essential').sum()):,} / "
        f"{len(genes[label]):,} essential overall" for label in MODELS
    )
    left.text(0, 4.1, counts_text, fontsize=10, va="center")
    names = substrates.substrate.drop_duplicates().tolist()
    ys = np.arange(len(names))
    for offset, label in ((-0.19, MODELS[0]), (0.19, MODELS[1])):
        frame = substrates[substrates.model.eq(label)].set_index("substrate").loc[names]
        right.barh(ys + offset, frame.growth.fillna(0.0), height=0.33, color=COLORS[label])
        right.scatter(np.zeros(len(frame)), ys + offset,
                      marker="o", s=14, c=[COLORS[label] if x == "zero_growth" else "none"
                                         for x in frame.call])
        for j, call in enumerate(frame.call):
            if call == "infeasible":
                right.scatter(0, ys[j] + offset, marker="x", s=24, c=COLORS[label],
                              clip_on=False)
    right.set_yticks(ys, names)
    right.invert_yaxis()
    right.set_xlabel("Optimal biomass flux (native objective of each model)")
    right.set_title("B  Single organic-carbon-source screen", loc="left", fontweight="bold")
    handles = [Patch(color=COLORS[label], label=label) for label in MODELS]
    right.legend(handles=handles, loc="lower right", frameon=False)
    right.text(1, -0.09, "● zero optimum    × infeasible",
               transform=right.transAxes, ha="right", fontsize=9)
    fig.suptitle("Bacillus subtilis 168  |  gapseq and AlphaGEM", fontsize=16, fontweight="bold")
    fig.savefig(out / "prediction_summary.png", dpi=180, bbox_inches="tight")
    fig.savefig(out / "prediction_summary.svg", bbox_inches="tight")
    plt.close(fig)


def run(paths, inputs, out):
    out.mkdir(parents=True, exist_ok=True)
    (out / "models").mkdir(exist_ok=True)
    media = pd.read_csv(inputs / "medium.csv", keep_default_na=False)
    menu = pd.read_csv(inputs / "substrates.csv")
    mapping = pd.read_csv(inputs / "uniprot_168.tsv", sep="\t", keep_default_na=False)
    proteome = mapping.set_index("Entry")["Gene Names"].to_dict()
    gene_frames, annotations = {}, {}
    substrates, summaries, controls, medium_sensitivity = [], [], [], []
    maintenance_sensitivity = pd.DataFrame()
    for label in MODELS:
        path = paths[label]
        reader = load_json_model if path.suffix == ".json" else read_sbml_model
        model = reader(str(path))
        model.solver = "glpk"
        exchanges = exchange_index(model)
        medium = read_medium(media, label, exchanges)
        controls.extend({"model": label, **row} for row in assay_controls(model, medium))
        for component in ("Molybdate", "Nickel"):
            exchange = media.loc[media.component.eq(component), label].iloc[0]
            glucose = menu.loc[menu.substrate.eq("D-Glucose")].iloc[0]
            base, _ = inputs_for_substrate(glucose, label, exchanges, medium)
            with model:
                model.medium = {rid: cap for rid, cap in base.items() if rid != exchange}
                removal_status, removal_growth = solve(model)
            medium_sensitivity.append({"model": label, "removed": component,
                                       "exchange": exchange, "status": removal_status,
                                       "growth": removal_growth})
        carbon = substrate_screen(model, label, menu, exchanges, medium)
        if label == "AlphaGEM" and "ATPM" in model.reactions:
            with model:
                model.reactions.ATPM.lower_bound = 0.0
                maintenance_sensitivity = substrate_screen(model, label, menu, exchanges, medium)
        glucose = menu.loc[menu.substrate.eq("D-Glucose")].iloc[0]
        model.medium, _ = inputs_for_substrate(glucose, label, exchanges, medium)
        status, baseline = solve(model)
        if status != "optimal" or baseline <= EPS:
            raise ValueError(f"{label} glucose positive control failed: {status}, {baseline}")
        model.id = f"Bsubtilis168_{label}_matched"
        write_sbml_model(model, str(out / "models" / f"{label}_matched.xml"))
        summaries.append(topology_summary(model, label, path, baseline))
        print(f"{label}: glucose growth {baseline:.7f}; screening {len(model.genes)} genes", flush=True)
        gene_frames[label] = gene_screen(model, label, baseline, proteome)
        annotations[label] = reaction_annotations(model)
        substrates.append(carbon)

    shared = gene_comparison(gene_frames)
    all_genes = pd.concat(gene_frames.values(), ignore_index=True)
    substrate_data = pd.concat(substrates, ignore_index=True)
    topology = pd.DataFrame(summaries)
    bigg = pd.DataFrame([
        {"bigg_id": value, **{label: "|".join(sorted(annotations[label].get(value, [])))
                              for label in MODELS}}
        for value in sorted(set(annotations[MODELS[0]]) | set(annotations[MODELS[1]]))
    ])
    for filename, data in {
        "gene_deletion": all_genes, "essential_genes": all_genes[all_genes.call.eq("essential")],
        "substrate_screen": substrate_data, "model_summary": topology,
        "controls": pd.DataFrame(controls), "bigg_reaction_overlap": bigg,
        "trace_medium_sensitivity": pd.DataFrame(medium_sensitivity),
        "alphagem_maintenance_sensitivity": maintenance_sensitivity,
        "shared_gene_comparison": shared.reset_index(),
    }.items():
        data.to_csv(out / f"{filename}.csv", index=False)
    plot(gene_frames, shared, substrate_data, out)
    (out / "parameters.json").write_text(json.dumps({
        "cobra_version": cobra.__version__, "solver": "glpk",
        "essential_fraction": ESSENTIAL_FRACTION, "growth_threshold": EPS,
        "carbon_flux_limit": CARBON_LIMIT,
    }, indent=2) + "\n")
    print(topology.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gapseq", type=Path, default=ROOT / "models/gapseq.xml")
    parser.add_argument("--alphagem", type=Path, default=ROOT / "models/alphagem.json")
    parser.add_argument("--inputs", type=Path, default=ROOT / "input")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    run(dict(zip(MODELS, (args.gapseq, args.alphagem), strict=True)), args.inputs, args.out)


if __name__ == "__main__":
    main()
