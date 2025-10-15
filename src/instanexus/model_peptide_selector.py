#!/usr/bin/env python

r"""Train and save a peptide selection model using Random Forest.
 _____  _______  _    _
|  __ \|__   __|| |  | |
| |  | |  | |   | |  | |
| |  | |  | |   | |  | |
| |__| |  | |   | |__| |
|_____/   |_|   |______|

__authors__ = Marco Reverenna
__copyright__ = Copyright 2025-2026
__research-group__ = DTU Biosustain (Multi-omics Network Analytics) and DTU Bioengineering
__date__ = 09 Oct 2025
__maintainer__ = Marco Reverenna
__email__ = marcor@dtu.dk
__status__ = Dev
"""


import re
import json
import pandas as pd
import numpy as np
import preprocessing as prep
from math import log2
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    precision_recall_curve,
    f1_score,
    average_precision_score,
)
import matplotlib.pyplot as plt
import joblib
import seaborn as sns

# side meaning where the protease cleaves
# residues meaning which amino acids it cleaves after (C-side) or before (N-side


def normalize_protease_list(p):
    """Normalize protease names (case-insensitive) to canonical lowercase aliases."""

    if p is None or (isinstance(p, float) and pd.isna(p)):
        return []

    if isinstance(p, str):
        toks = re.split(r"[;,\|\s]+", p.strip().lower())
    elif isinstance(p, (list, tuple, set)):
        toks = [str(x).lower() for x in p]
    else:
        toks = [str(p).lower()]

    alias = {
        "trypsin": "trypsin",
        "lysc": "lysc",
        "gluc": "gluc",
        "chymotrypsin": "chymotrypsin",
        "chymo": "chymotrypsin",
        "chemo": "chymotrypsin",
        "elastase": "elastase",
        "papain": "papain",
        "legumain": "legumain",
        "protk": "protk",
        "proteinasek": "protk",
        "krakatoa": "krakatoa",
        "thermo": "thermo",
        "thermolysin": "thermo",
        "vesuvius": "vesuvius",
    }

    return [alias.get(t, t) for t in toks if t]


def load_protease_rules(json_path):
    """Load protease cleavage rules from a JSON file and convert lists to sets."""
    with open(json_path, "r") as f:
        rules = json.load(f)
    for rule in rules.values():
        rule["residues"] = set(rule["residues"])
        rule["no_cut_if_next_is"] = set(rule["no_cut_if_next_is"])
        rule["no_cut_if_prev_is"] = set(rule["no_cut_if_prev_is"])
    return rules


def cterm_matches_any(seq, prots, protease_rules):
    """Check if peptide C-terminal matches any protease cleavage rule."""
    if not seq:
        return 0

    last = seq[-1]  # cleavage site residue
    prev_res = seq[-2] if len(seq) > 1 else None  # residue before cleavage

    for p in prots:
        rule = protease_rules.get(p)
        if not rule or rule["side"] != "C":
            continue

        if last in rule["residues"]:
            # inhibit if previous residue blocks cleavage
            if prev_res and prev_res in rule["no_cut_if_prev_is"]:
                continue

            # inhibit if proline follows cleavage site
            if prev_res == "P" and "P" in rule["no_cut_if_next_is"]:
                continue

            return 1

    return 0


def nterm_matches_any(seq, prots, protease_rules):
    """Check if peptide N-terminal matches any protease cleavage rule."""
    if not seq:
        return 0

    first = seq[0]  # residue after cleavage
    next_res = seq[1] if len(seq) > 1 else None

    for p in prots:
        rule = protease_rules.get(p)
        if not rule or rule["side"] != "N":
            continue

        # cleavage occurs before this residue
        if first in rule["residues"]:
            # block if the previous residue inhibits cleavage
            # (conceptually: residue before cleavage, i.e. at the C-term of previous peptide)
            if "P" in rule["no_cut_if_prev_is"]:
                continue

            # block if the next residue inhibits cleavage
            if next_res and next_res in rule["no_cut_if_next_is"]:
                continue

            return 1

    return 0


def internal_expected_sites_min(seq, prots, protease_rules):
    """Count minimal expected internal cleavage sites across proteases,
    respecting both no_cut_if_prev_is and no_cut_if_next_is rules."""
    if not seq or len(seq) < 2 or not prots:
        return 0

    counts = []

    for p in prots:
        rule = protease_rules.get(p)
        if not rule:
            continue

        cnt = 0

        # C-terminal cleavage protease
        if rule["side"] == "C":
            for i, a in enumerate(seq[:-1]):  # cleavage after position i
                next_res = seq[i + 1] if i + 1 < len(seq) else None
                prev_res = seq[i - 1] if i > 0 else None

                # valid cleavage target?
                if a not in rule["residues"]:
                    continue

                # block if previous residue inhibits cleavage
                if prev_res and prev_res in rule["no_cut_if_prev_is"]:
                    continue

                # block if next residue inhibits cleavage
                if next_res and next_res in rule["no_cut_if_next_is"]:
                    continue

                cnt += 1

        # N-terminal cleavage protease
        elif rule["side"] == "N":
            for i, a in enumerate(seq[1:], start=1):  # cleavage before position i
                next_res = seq[i + 1] if i + 1 < len(seq) else None
                prev_res = seq[i - 1] if i > 0 else None

                # valid cleavage target?
                if a not in rule["residues"]:
                    continue

                # block if previous residue inhibits cleavage
                if prev_res and prev_res in rule["no_cut_if_prev_is"]:
                    continue

                # block if next residue inhibits cleavage
                if next_res and next_res in rule["no_cut_if_next_is"]:
                    continue

                cnt += 1

        counts.append(cnt)

    return min(counts) if counts else 0


def proline_block_at_cterm(seq):
    """Return 1 if residue before C-term is proline, else 0."""
    if not seq:
        return 0
    return int(len(seq) >= 2 and seq[-2] == "P")


def max_repeat_length(seq):
    """Return length of the longest consecutive amino acid repeat."""
    if not isinstance(seq, str) or len(seq) == 0:
        return 0
    matches = list(re.finditer(r"(.)\1+", seq))
    return max((len(m.group(0)) for m in matches), default=1)


def repeat_ratio(seq):
    """Compute fraction of repeated amino acids over sequence length."""
    if not isinstance(seq, str) or len(seq) == 0:
        return 0
    total_repeats = sum(len(m.group(0)) for m in re.finditer(r"(.)\1+", seq))
    return total_repeats / len(seq)


def seq_entropy(seq):
    """Compute Shannon entropy of amino acid composition."""
    if not isinstance(seq, str) or len(seq) == 0:
        return 0
    freqs = [seq.count(a) / len(seq) for a in set(seq)]
    return -sum(p * log2(p) for p in freqs)


def load_aa_properties(json_path):
    """Load amino acid properties (hydrophobicity, mass, etc.) from JSON."""
    with open(json_path, "r") as f:
        return json.load(f)


def peptide_props(seq, aa_properties):
    """Calculate hydrophobicity, mass stats, and basic residue fraction."""
    if not seq or not isinstance(seq, str) or len(seq) == 0:
        return pd.Series(
            {"mean_hydro": 0, "mean_mass": 0, "mass_std": 0, "frac_basic": 0}
        )

    vals_h = [aa_properties.get(a, {"hydro": 0})["hydro"] for a in seq]
    vals_m = [aa_properties.get(a, {"mass": 0})["mass"] for a in seq]
    frac_basic = sum(seq.count(a) for a in "KRH") / len(seq)

    return pd.Series(
        {
            "mean_hydro": np.mean(vals_h),
            "mean_mass": np.mean(vals_m),
            "mass_std": np.std(vals_m),
            "frac_basic": frac_basic,
        }
    )


def build_reference_free_features(df, aa_properties, protease_rules):
    """Generate all peptide-level numerical features for model training."""

    df = df.copy()
    df["seq_length"] = df["cleaned_preds"].str.len()
    df["has_special"] = (
        df["cleaned_preds"].str.contains(r"[^A-Z]", regex=True).astype(int)
    )
    df["first_aa"] = df["cleaned_preds"].str[0].astype("category").cat.codes
    df["last_aa"] = df["cleaned_preds"].str[-1].astype("category").cat.codes

    # Physicochemical properties
    df = pd.concat(
        [df, df["cleaned_preds"].apply(lambda s: peptide_props(s, aa_properties))],
        axis=1,
    )

    df["has_proline"] = df["cleaned_preds"].str.contains("P").astype(int)
    df["has_cysteine"] = df["cleaned_preds"].str.contains("C").astype(int)
    df["contains_nq"] = df["cleaned_preds"].str.contains("[NQ]").astype(int)
    df["has_repeated_aa"] = df["cleaned_preds"].str.contains(r"(.)\1").astype(int)
    df["max_repeat_len"] = df["cleaned_preds"].apply(max_repeat_length)
    df["repeat_ratio"] = df["cleaned_preds"].apply(repeat_ratio)
    df["seq_entropy"] = df["cleaned_preds"].apply(seq_entropy)

    prots_list = df["protease"].apply(normalize_protease_list)

    df["cterm_matches_protease"] = [
        cterm_matches_any(s, p, protease_rules)
        for s, p in zip(df["cleaned_preds"].fillna(""), prots_list)
    ]
    df["nterm_matches_protease"] = [
        nterm_matches_any(s, p, protease_rules)
        for s, p in zip(df["cleaned_preds"].fillna(""), prots_list)
    ]
    df["internal_expected_sites_min"] = [
        internal_expected_sites_min(s, p, protease_rules)
        for s, p in zip(df["cleaned_preds"].fillna(""), prots_list)
    ]

    df["proline_block_at_cterm"] = (
        df["cleaned_preds"].fillna("").apply(proline_block_at_cterm)
    )
    df["protease"] = df["protease"].astype("category").cat.codes

    return df


def train_model(df, reference_seq, model_path, aa_properties, protease_rules):
    """Train Random Forest classifier and save model with optimal threshold."""
    df["mapped"] = df["cleaned_preds"].apply(
        lambda x: int(isinstance(x, str) and x in reference_seq)
    )
    df = build_reference_free_features(df, aa_properties, protease_rules)

    exclude = ["experiment_name", "scan_number", "preds", "cleaned_preds"]
    feature_cols = [c for c in df.columns if c not in exclude and c != "mapped"]

    x = df[feature_cols]
    y = df["mapped"].astype(int)
    x_train, x_test, y_train, y_test = train_test_split(
        x, y, test_size=0.3, stratify=y, random_state=42
    )

    model = RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1)
    model.fit(x_train, y_train)
    y_scores = model.predict_proba(x_test)[:, 1]

    precision, recall, thresholds = precision_recall_curve(y_test, y_scores)
    f1_scores = [f1_score(y_test, (y_scores >= t).astype(int)) for t in thresholds]
    best_idx = max(range(len(f1_scores)), key=f1_scores.__getitem__)
    best_threshold = thresholds[best_idx]
    ap = average_precision_score(y_test, y_scores)

    joblib.dump(
        {"model": model, "threshold": best_threshold, "features": feature_cols},
        model_path,
    )

    print(f"\n Model trained and saved to: {model_path}")
    print(f"Best F1 threshold: {best_threshold:.3f}")
    print(f"Average precision (AP): {ap:.3f}")

    # Return metrics for external plotting
    return {
        "precision": precision,
        "recall": recall,
        "f1_scores": f1_scores,
        "thresholds": thresholds,
        "best_idx": best_idx,
        "ap": ap,
    }


def plot_precision_recall(metrics, output_dir, filename="precision_recall_curve.svg"):
    """Seaborn Precision–Recall curve."""
    sns.set_theme(style="white", font_scale=1)
    plt.figure(figsize=(7, 5))

    recall = metrics["recall"]
    precision = metrics["precision"]
    best_idx = metrics["best_idx"]
    ap = metrics["ap"]

    sns.lineplot(
        x=recall, y=precision, color="#2E86AB", linewidth=1, label=f"AP = {ap:.2f}"
    )
    plt.scatter(
        recall[best_idx],
        precision[best_idx],
        color="#D64550",
        s=80,
        label="Max F1",
        zorder=5,
    )

    plt.xlabel("Recall", fontsize=12)
    plt.ylabel("Precision", fontsize=12)
    plt.title("Precision–Recall Curve", fontsize=13, pad=15)
    plt.legend(frameon=False)
    sns.despine()
    plt.tight_layout()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / filename
    plt.savefig(out_path, format="svg", dpi=300)
    plt.close()
    print(f"Saved: {out_path}")


def get_sample_metadata(run, chain, json_path):
    with open(json_path, "r") as f:
        all_meta = json.load(f)

    if run not in all_meta:
        raise ValueError(f"Run '{run}' not found in metadata.")

    entries = all_meta[run]

    for entry in entries:
        if entry["chain"] == chain:
            return entry

    raise ValueError(f"No metadata found for run '{run}' with chain '{chain}'.")


def main():
    """Main execution for training and saving the peptide selection model."""
    try:
        # works if you are in a script: __file__ exists
        BASE_DIR = Path(__file__).resolve().parents[1]
    except NameError:
        # works if you are in a notebook: __file__ does not exist
        BASE_DIR = Path().resolve()
        # go up until the project folder
        while BASE_DIR.name != "InstaNexus" and BASE_DIR != BASE_DIR.parent:
            BASE_DIR = BASE_DIR.parent

    JSON_DIR = BASE_DIR / "json"
    INPUT_DIR = BASE_DIR / "inputs"
    FASTA_DIR = BASE_DIR / "fasta"
    FIGURE_DIR = BASE_DIR / "figures"

    run = "bsa"
    chain = ""

    meta = get_sample_metadata(run, chain, json_path=JSON_DIR / "sample_metadata.json")

    protein = meta["protein"]
    chain = meta["chain"]
    proteases = meta["proteases"]

    aa_props = load_aa_properties(JSON_DIR / "aa_properties.json")
    protease_rules = load_protease_rules(JSON_DIR / "protease_rules.json")

    protein_norm = prep.normalize_sequence(protein)
    df = pd.read_csv(INPUT_DIR / f"{run}.csv")

    df["protease"] = df["experiment_name"].apply(
        lambda name: prep.extract_protease(name, proteases)
    )

    df = prep.clean_dataframe(df)

    df["cleaned_preds"] = df["preds"].apply(prep.remove_modifications)

    cleaned_psms = df["cleaned_preds"].tolist()

    filtered_psms = prep.filter_contaminants(
        cleaned_psms, run, FASTA_DIR / "contaminants.fasta"
    )

    df = df[df["cleaned_preds"].isin(filtered_psms)]

    df["mapped"] = df["cleaned_preds"].apply(
        lambda x: int(isinstance(x, str) and x in protein_norm)
    )

    model_path = BASE_DIR / "peptide_selector.pkl"
    metrics = train_model(df, protein, model_path, aa_props, protease_rules)
    plot_precision_recall(metrics, FIGURE_DIR)


if __name__ == "__main__":
    main()
