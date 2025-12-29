# logic.py
import requests
from collections import defaultdict
from typing import Dict

# ------------------------------
# IMPORT ONLY FUNCTIONS YOU USE
# ------------------------------
from infinity_stone_2 import (
    gene_symbol_to_uniprot,
    get_e2s_from_reactome,
    get_e2s_from_biogrid,
    get_go_terms_from_uniprot,
    score_e2_type
)

from infinity_stone_4 import get_accessible_lysines
from infinity_stone_5 import analyze_substrate
from infinity_stone_7 import infer_local_dub_activity
from infinity_stone_8 import biochemical_linkage_override
from infinity_stone_9 import apply_temporal_gate_explicit
from infinity_stone_10 import (
    find_prior_files,
    load_priors,
    map_linkage_to_process
)

from infinity_stone_1 import (
    get_ubibrowser_term_from_uniprot,
    get_e3_class_from_ubibrowser,
    E3_POLYUB_SCORE
)

# ------------------------------
# HELPERS
# ------------------------------
def canonicalize_linkage(l):
    if l in {"M1", "linear"}:
        return "M1_linear"
    return l

def resolve_e3_uniprot_and_gene(e3_input):
    if e3_input.upper().startswith(("Q", "P")):
        return e3_input, None
    return gene_symbol_to_uniprot(e3_input), e3_input

# ------------------------------
# MAIN PIPELINE (PURE FUNCTION)
# ------------------------------
def run_pipeline(
    substrate_uniprot: str,
    e3_input: str,
    cell_state: str,
    cell_type: str,
    top_n_e2: int = 3
) -> Dict:

    # ---- Resolve E3 ----
    e3_uniprot, e3_gene = resolve_e3_uniprot_and_gene(e3_input)

    poly_prior = 0.3
    if e3_uniprot:
        term = get_ubibrowser_term_from_uniprot(e3_uniprot)
        cls = get_e3_class_from_ubibrowser(term) if term else None
        poly_prior = E3_POLYUB_SCORE.get(cls, poly_prior)

    # ---- E2s ----
    e2s = get_e2s_from_reactome() or get_e2s_from_biogrid(e3_gene)
    top_e2s = sorted(list(e2s))[:top_n_e2]

    # ---- Substrate ----
    lys_cat, lys_det = get_accessible_lysines(substrate_uniprot)
    residence = analyze_substrate(substrate_uniprot, e3_gene)

    # ---- DUB ----
    local_dub_state, local_dub_evidence = infer_local_dub_activity(
        substrate_uniprot, e3_uniprot
    )

    # ---- Biochemical aggregation ----
    agg = defaultdict(float)

    for e2 in top_e2s:
        res = biochemical_linkage_override(
            substrate_uniprot,
            e3_gene or "",
            e2_gene=e2,
            local_dubs=list(local_dub_evidence.get("DUBs", {}).keys())
        )
        for k, v in res.get("scores", {}).items():
            agg[canonicalize_linkage(k)] += max(v, 0)

    if not agg:
        agg = {"K48": 0.5, "K63": 0.3, "MONO": 0.2}

    total = sum(agg.values())
    probs = {k: v / total for k, v in agg.items()}

    # ---- Temporal gate ----
    gated, _ = apply_temporal_gate_explicit(probs, cell_state)
    final = gated if gated else probs

    dominant_linkage = max(final, key=final.get)

    # ---- Process inference ----
    priors = load_priors(find_prior_files())
    process_df, _ = map_linkage_to_process(
        cell_type, dominant_linkage, priors
    )

    return {
        "dominant_linkage": dominant_linkage,
        "linkage_probs": final,
        "processes": process_df.to_dict(orient="records"),
        "lysine_category": lys_cat,
        "local_dub_state": local_dub_state
    }
