"""
INFINITY STONE 10
=================
Cell-type + ubiquitin-linkage → process bias mapper (BN-friendly)
"""

import os
import glob
import pandas as pd
from typing import Dict, Tuple

# =========================================================
# CANONICAL USER CHOICES
# =========================================================

CELL_TYPE_CHOICES = [
    "Neuronal cells",
    "Glial cells",
    "Endocrine cells",
    "Squamous epithelial cells",
    "Pigment cells",
    "Ciliated cells",
    "Specialized epithelial cells",
    "Glandular epithelial cells",
    "Germ cells",
    "Trophoblast cells",
    "Muscle cells",
    "Endothelial and mural cells",
    "Mesenchymal cells",
    "Blood and immune cells",
    "Stem and proliferating cells",
]

UBIQUITIN_LINKAGES = [
    "MONO",
    "K48", "K63", "K11", "K29",
    "K6", "K27", "K33",
    "M1",
]

# =========================================================
# LINKAGE → PROCESS LIKELIHOODS  P(linkage | process)
# =========================================================

LINKAGE_LIKELIHOOD = {
    "MONO": {"proteasome": 0.01, "lysosome": 0.35, "autophagy": 0.15, "recycling": 0.30, "escrt": 0.19},

    "K48": {"proteasome": 0.75, "lysosome": 0.05, "autophagy": 0.05, "recycling": 0.02, "escrt": 0.02},
    "K11": {"proteasome": 0.45, "lysosome": 0.05, "autophagy": 0.05, "recycling": 0.02, "escrt": 0.02},
    "K63": {"proteasome": 0.05, "lysosome": 0.35, "autophagy": 0.35, "recycling": 0.10, "escrt": 0.10},
    "K29": {"proteasome": 0.05, "lysosome": 0.35, "autophagy": 0.20, "recycling": 0.15, "escrt": 0.10},
    "K6":  {"proteasome": 0.10, "lysosome": 0.20, "autophagy": 0.10, "recycling": 0.10, "escrt": 0.05},
    "K27": {"proteasome": 0.05, "lysosome": 0.20, "autophagy": 0.15, "recycling": 0.10, "escrt": 0.05},
    "K33": {"proteasome": 0.05, "lysosome": 0.10, "autophagy": 0.10, "recycling": 0.15, "escrt": 0.05},
    "M1":  {"proteasome": 0.02, "lysosome": 0.02, "autophagy": 0.02, "recycling": 0.02, "escrt": 0.02},
}

# Normalize likelihood vectors (unchanged behavior)
for lk, vec in LINKAGE_LIKELIHOOD.items():
    s = sum(vec.values())
    for p in vec:
        vec[p] /= s

# =========================================================
# LINKAGE EDITABILITY (DUB sensitivity)
# =========================================================

LINKAGE_EDITABILITY = {
    "MONO": 0.15,
    "K48": 0.8,
    "K11": 0.6,
    "K63": 0.6,
    "K29": 0.5,
    "K6":  0.4,
    "K27": 0.4,
    "K33": 0.3,
    "M1":  0.2,
}

EPSILON_PRIOR = 1e-6

# =========================================================
# UTILITIES
# =========================================================

def find_prior_files(pattern="*_cell_type_priors.csv") -> Dict[str, str]:
    files = glob.glob(pattern)
    mapping = {}
    for f in files:
        process = os.path.basename(f).split("_cell_type_priors")[0].lower()
        mapping[process] = f
    return mapping


def load_priors(prior_files: Dict[str, str]) -> Dict[str, pd.Series]:
    priors = {}
    for process, path in prior_files.items():
        df = pd.read_csv(path)
        s = pd.Series(df["Prior"].values, index=df["Cell_Type"].astype(str))
        priors[process] = s.astype(float).clip(0.0, 1.0)
    return priors

# =========================================================
# CORE MAPPING FUNCTION
# =========================================================

def map_linkage_to_process(
    cell_type: str,
    linkage: str,
    priors_dict: Dict[str, pd.Series],
    dub_process_name: str = "deubiquitination",
) -> Tuple[pd.DataFrame, Dict]:

    linkage = linkage.upper()
    processes = sorted(priors_dict.keys())

    # P(process | cell)
    p_cell = {
        p: float(priors_dict[p].get(cell_type, EPSILON_PRIOR))
        for p in processes
    }

    # P(linkage | process)
    p_l_given_proc = {
        p: LINKAGE_LIKELIHOOD.get(linkage, {}).get(p, EPSILON_PRIOR)
        for p in processes
    }

    dub_series = priors_dict.get(dub_process_name, pd.Series(dtype=float))
    dub_prior = float(dub_series.get(cell_type, 0.0))
    editability = LINKAGE_EDITABILITY.get(linkage, 0.4)

    def is_degradative(proc):
        return any(k in proc.lower() for k in ("proteasome", "lysosome", "autophagy"))

    rows = []
    for p in processes:
        survival = (
            1.0 - dub_prior * editability
            if is_degradative(p)
            else 1.0 + 0.5 * dub_prior * editability
        )
        survival = max(0.0, survival)

        unnorm = p_cell[p] * p_l_given_proc[p] * survival

        rows.append({
            "Process": p,
            "Prior_given_cell": p_cell[p],
            "P_linkage_given_process": p_l_given_proc[p],
            "Survival_factor": round(survival, 4),
            "Unnormalized_score": unnorm,
        })

    df = pd.DataFrame(rows)

    if df.empty or df["Unnormalized_score"].sum() == 0:
        df["Posterior"] = [
            p_cell[p] / sum(p_cell.values())
            for p in df["Process"]
        ]
    else:
        total = df["Unnormalized_score"].sum()
        df["Posterior"] = df["Unnormalized_score"] / total

    df = df.sort_values("Posterior", ascending=False).reset_index(drop=True)

    metadata = {
        "cell_type": cell_type,
        "linkage": linkage,
        "dub_prior": dub_prior,
        "editability": editability,
    }

    return (
        df[
            [
                "Process",
                "Posterior",
                "Prior_given_cell",
                "P_linkage_given_process",
                "Survival_factor",
            ]
        ],
        metadata,
    )
