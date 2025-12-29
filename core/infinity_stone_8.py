# =========================
# E2-driven linkage specificity
# =========================
E2_LINKAGE_RULES = {
    "UBE2S":  {"K11": 0.95},
    "UBE2N":  {"K63": 0.9},
    "UBE2V1": {"K63": 0.9},
    "UBE2V2": {"K63": 0.9},
    "UBE2D":  {"K48": 0.5, "K63": 0.5},
}

# =========================
# E3-driven linkage rules
# =========================
E3_LINKAGE_RULES = {
    "APC/C":  {"K11": 1.0},
    "ANAPC":  {"K11": 1.0},
    "BRCA1":  {"K6": 1.0},
    "BARD1":  {"K6": 1.0},
    "HOIP":   {"M1": 1.0},
    "RNF31":  {"M1": 1.0},
    "TRIM13": {"K11": 0.8},
    "TRIM25": {"K63": 0.7, "K48": 0.3},
}

# =========================
# DUB anti-linkage pressure
# =========================
DUB_ANTI_LINKAGE = {
    "OTULIN": {"M1": -1.0},
    "CEZANNE": {"K11": -0.8},
    "TRABID": {"K29": -0.8, "K33": -0.8},
    "CYLD": {"K63": -0.6},
}

# =========================
# Helpers
# =========================
def normalize_gene(g):
    return g.upper().replace("-", "").replace("_", "")

def resolve_e3_family(e3_gene):
    g = normalize_gene(e3_gene)
    if g.startswith("ANAPC"):
        return "APC/C"
    if g in {"BRCA1", "BARD1"}:
        return g
    if g in {"RNF31", "HOIP"}:
        return "HOIP"
    return g

# =========================
# MAIN ENTRY POINT
# =========================
def biochemical_linkage_override(
    substrate_ac,
    e3_gene,
    e2_gene=None,
    local_dubs=None
):
    linkage_scores = {}

    e3_key = resolve_e3_family(e3_gene)
    if e3_key in E3_LINKAGE_RULES:
        for link, score in E3_LINKAGE_RULES[e3_key].items():
            linkage_scores[link] = linkage_scores.get(link, 0) + score

    if e2_gene:
        e2_key = normalize_gene(e2_gene)
        if e2_key in E2_LINKAGE_RULES:
            for link, score in E2_LINKAGE_RULES[e2_key].items():
                linkage_scores[link] = linkage_scores.get(link, 0) + score

    if local_dubs:
        for dub in local_dubs:
            dub_key = normalize_gene(dub)
            if dub_key in DUB_ANTI_LINKAGE:
                for link, penalty in DUB_ANTI_LINKAGE[dub_key].items():
                    linkage_scores[link] = linkage_scores.get(link, 0) + penalty

    if not linkage_scores:
        return {
            "mode": "system_level_only",
            "allowed_linkages": ["K48", "K63", "mixed"],
            "scores": {}
        }

    linkage_scores = {
        k: round(max(min(v, 1.0), -1.0), 3)
        for k, v in linkage_scores.items()
    }

    return {
        "mode": "biochemical_override",
        "dominant_linkage": max(linkage_scores, key=linkage_scores.get),
        "scores": linkage_scores
    }
