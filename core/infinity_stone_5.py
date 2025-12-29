import requests
import re
from bs4 import BeautifulSoup

# ==================================================
# CONSTANTS
# ==================================================

# ---- DegPred ----
DEGPRED_URL = "http://degron.phasep.pro/detail/{}"

# ---- BioGRID ----
BIOGRID_URL = "https://webservice.thebiogrid.org/interactions"
BIOGRID_ACCESS_KEY = "1ab875d99061e9ae48d022f346846e3b"
BIOGRID_SPECIES = 9606

# ---- Reactome ----
REACTOME_BASE = "https://reactome.org/ContentService"
REACTOME_SPECIES = "Homo sapiens"

# ==================================================
# GO CONTEXT MAP (RESTORED)
# ==================================================

GO_CONTEXT_MAP = {

    # =========================
    # ER-ASSOCIATED DEGRADATION
    # =========================
    "erad": {
        "GO:0036503",
        "GO:0097466",
        "GO:1904263",
        "GO:0045047",
        "GO:0016567",
        "GO:0043161",
        "GO:0005783",
        "GO:0031965",
    },

    # =========================
    # EGAD (ENDOSOME / GOLGI)
    # =========================
    "egad": {
        "GO:0005769",
        "GO:0005770",
        "GO:0032585",
        "GO:0005794",
        "GO:0006511",
        "GO:0043161",
        "GO:0016567",
        "GO:0006897",
        "GO:0006898",
        "GO:0015031",
    },

    # =========================
    # SIGNALING
    # =========================
    "signaling": {
        "GO:0007165",
        "GO:0023052",
        "GO:0006955",
        "GO:0007249",
        "GO:0043123",
        "GO:0000165",
        "GO:0070371",
        "GO:0007169",
        "GO:0007186",
        "GO:0018105",
        "GO:0006468",
    }
}

# ==================================================
# NEGATIVE EXCLUSION RULES (RESTORED)
# ==================================================

GO_NEGATIVE_EXCLUSION = {

    "erad": {
        "GO:0005769",
        "GO:0005770",
        "GO:0032585",
        "GO:0005794",
        "GO:0005886",
        "GO:0009897",
    },

    "egad": {
        "GO:0036503",
        "GO:0097466",
        "GO:0045047",
        "GO:0005783",
        "GO:0031965",
    },

    "signaling": {
        "GO:0006511",
        "GO:0043161",
        "GO:0005198",
        "GO:0005737",
    }
}

CONTEXT_PRIORITY = ["erad", "egad", "signaling"]

# ==================================================
# BIOGRID INTERACTION COUNT
# ==================================================

def get_biogrid_interaction_count(protein):
    params = {
        "geneList": protein,
        "format": "json",
        "accessKey": BIOGRID_ACCESS_KEY,
        "organism": BIOGRID_SPECIES,
        "includeInteractors": "true",
        "interSpeciesExcluded": "true"
    }

    try:
        r = requests.get(BIOGRID_URL, params=params, timeout=30)
        r.raise_for_status()
        return len(r.json())
    except Exception:
        return 0

# ==================================================
# REACTOME → GO
# ==================================================

def reactome_pathway_to_go_accession(pathway_id):
    url = f"{REACTOME_BASE}/data/query/{pathway_id}"
    r = requests.get(url, headers={"accept": "application/json"}, timeout=20)
    r.raise_for_status()

    go_bp = r.json().get("goBiologicalProcess")
    if not go_bp or not go_bp.get("accession"):
        return None

    return f"GO:{go_bp['accession']}"

def get_reactome_context_via_reactions(uniprot_id):
    url = (
        f"{REACTOME_BASE}/data/mapping/UniProt/"
        f"{uniprot_id}/reactions?species=9606"
    )

    r = requests.get(url, headers={"accept": "application/json"}, timeout=20)
    r.raise_for_status()
    reactions = r.json()

    pathway_ids = set()
    for rxn in reactions:
        if rxn.get("dbId"):
            pathway_ids.add(rxn["dbId"])

    go_terms = set()
    for pid in pathway_ids:
        go_id = reactome_pathway_to_go_accession(pid)
        if go_id:
            go_terms.add(go_id)

    context_hits = assign_context_with_exclusions(go_terms)
    resolved_context, confidence = resolve_context(context_hits)

    return {
        "reaction_count": len(reactions),
        "pathway_count": len(pathway_ids),
        "go_terms": sorted(go_terms),
        "context_hits": context_hits,
        "pathway_context": resolved_context,
        "context_confidence": confidence
    }

# ==================================================
# GO CONTEXT RESOLUTION
# ==================================================

def assign_context_with_exclusions(go_terms):
    context_hits = {}

    for ctx, positives in GO_CONTEXT_MAP.items():
        pos_hits = go_terms & positives
        neg_hits = go_terms & GO_NEGATIVE_EXCLUSION.get(ctx, set())

        if pos_hits:
            context_hits[ctx] = {
                "positive_hits": sorted(pos_hits),
                "negative_hits": sorted(neg_hits),
                "score": len(pos_hits) - len(neg_hits)
            }

    return context_hits

def resolve_context(context_hits):
    if not context_hits:
        return "other", "low"

    for ctx in CONTEXT_PRIORITY:
        if ctx in context_hits:
            if not context_hits[ctx]["negative_hits"]:
                return ctx, "high"
            return ctx, "conflicted"

    return "other", "low"

# ==================================================
# DEGPRED PARSER (UNCHANGED)
# ==================================================

def parse_degpred(uniprot_id):
    url = DEGPRED_URL.format(uniprot_id)
    r = requests.get(url, timeout=30)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")

    result = {
        "degpred_regions": [],
        "elm_regions": [],
        "disorder_regions": []
    }

    for script in soup.find_all("script"):
        if not script.string:
            continue
        text = script.string

        if "var degpredregions" in text:
            matches = re.findall(
                r"\{x:\s*\[(\d+),.*?\]\[i\],\s*y:\s*\[(\d+),.*?\]\[i\],.*?\['(.*?)'",
                text
            )
            for m in matches:
                result["degpred_regions"].append({
                    "start": int(m[0]),
                    "end": int(m[1]),
                    "sequence": m[2]
                })

        if "var elmregions" in text:
            matches = re.findall(
                r"\{x:\s*\[(\d+),.*?\]\[i\],\s*y:\s*\[(\d+),.*?\]\[i\],.*?\['(.*?)'\]",
                text
            )
            for m in matches:
                result["elm_regions"].append({
                    "start": int(m[0]),
                    "end": int(m[1]),
                    "motif": m[2]
                })

        if "var espritzregion" in text:
            matches = re.findall(
                r"\{x:\s*\[(\d+)\]\[i\],\s*y:\s*\[(\d+)\]\[i\]",
                text
            )
            for m in matches:
                result["disorder_regions"].append({
                    "start": int(m[0]),
                    "end": int(m[1])
                })

    return result

# ==================================================
# RESIDENCE TIME INFERENCE (UNCHANGED)
# ==================================================

def overlaps(region, disorder_regions):
    return any(
        region["start"] <= d["end"] and region["end"] >= d["start"]
        for d in disorder_regions
    )

def analyze_residence_time(parsed, biogrid_count, reactome_ctx):
    predicted_degrons = len(parsed["degpred_regions"])
    degron_diversity = len(set(e.get("motif") for e in parsed["elm_regions"]))
    degrons_in_disorder = any(
        overlaps(d, parsed["disorder_regions"])
        for d in parsed["degpred_regions"]
    )

    score = 0
    if predicted_degrons >= 2:
        score += 1
    if degron_diversity >= 2:
        score += 1
    if degrons_in_disorder:
        score += 1
    if biogrid_count >= 10:
        score += 1
    if (
        reactome_ctx["pathway_context"] in {"egad", "erad"}
        and reactome_ctx["context_confidence"] == "high"
    ):
        score += 1

    return {
        "residence_score": score,
        "Substrate_residence_time": "long" if score >= 3 else "short"
    }

# ==================================================
# MASTER ENTRY POINT
# ==================================================

def analyze_substrate(uniprot_id, gene_symbol):
    parsed = parse_degpred(uniprot_id)
    biogrid_count = get_biogrid_interaction_count(gene_symbol)
    reactome_ctx = get_reactome_context_via_reactions(uniprot_id)

    residence = analyze_residence_time(parsed, biogrid_count, reactome_ctx)

    return {
        "uniprot_id": uniprot_id,
        "gene_symbol": gene_symbol,
        "degpred": parsed,
        "biogrid_interactions": biogrid_count,
        "reactome_context": reactome_ctx,
        "residence_time_analysis": residence
    }
