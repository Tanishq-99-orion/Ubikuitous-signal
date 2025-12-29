import requests
import re
import json
from bs4 import BeautifulSoup

# =========================
# DUB FAMILY BIAS (UNCHANGED)
# =========================
DUB_FAMILY_BIAS = {
    "USP": {"bias": "stabilizing", "score": 0.7},
    "OTU": {"bias": "chain_specific", "score": 0.6},
    "JAMM": {"bias": "endosomal", "score": 0.8},
    "MINDY": {"bias": "k48_trimming", "score": 0.85},
    "UCH": {"bias": "recycling", "score": 0.4},
    "Josephin": {"bias": "quality_control", "score": 0.6},
    "Other": {"bias": "unknown", "score": 0.3},
}

# =========================
# SUBCELLULAR WEIGHTS (HPA)
# =========================
SUBCELLULAR_BASE_WEIGHTS = {
    "nucleoplasm": 1.0,
    "nuclear bodies": 0.9,
    "nuclear speckles": 0.85,
    "nucleoli": 0.8,
    "endosomes": 1.0,
    "lysosomes": 0.9,
    "golgi apparatus": 0.7,
    "endoplasmic reticulum": 0.6,
    "vesicles": 0.6,
    "cytosol": 0.6,
    "centrosome": 0.5,
    "mitochondria": 0.4,
}

LOCATION_GROUP_MAP = {
    "nucleoplasm": "nucleus",
    "nuclear bodies": "nucleus",
    "nuclear speckles": "nucleus",
    "nucleoli": "nucleus",
    "endosomes": "endosome",
    "lysosomes": "endosome",
    "golgi apparatus": "golgi",
    "endoplasmic reticulum": "ERAD",
    "vesicles": "endosome",
    "cytosol": "cytosol",
    "centrosome": "cytosol",
    "mitochondria": "mitochondria",
}

SUBCELLULAR_DUB_PRIORS = {
    "nucleus": 0.8,
    "endosome": 0.85,
    "golgi": 0.6,
    "cytosol": 0.5,
    "ERAD": 0.25,
    "mitochondria": 0.3,
    "unknown": 0.3,
}

# =========================
# UniProt helpers
# =========================
def get_canonical_gene_symbol(uniprot_ac):
    try:
        r = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json",
            timeout=15
        )
        r.raise_for_status()
        for gene in r.json().get("genes", []):
            if gene.get("geneName", {}).get("value"):
                return gene["geneName"]["value"]
    except Exception:
        pass
    return None


def gene_symbol_to_ensembl(symbol):
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}"
    try:
        r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=15)
        r.raise_for_status()
        for entry in r.json():
            if entry.get("type") == "gene" and entry.get("id", "").startswith("ENSG"):
                return entry["id"]
    except Exception:
        pass
    return None

# =========================
# HPA subcellular inference
# =========================
def infer_subcellular_location_hpa_from_symbol(gene_symbol):
    ensg = gene_symbol_to_ensembl(gene_symbol)
    if not ensg:
        return "unknown", {}, {}

    try:
        r = requests.get(
            f"https://www.proteinatlas.org/{ensg}/subcellular",
            timeout=20
        )
        r.raise_for_status()
        text = BeautifulSoup(r.text, "html.parser").get_text(" ").lower()
    except Exception:
        return "unknown", {}, {}

    detected = set()
    for loc in SUBCELLULAR_BASE_WEIGHTS:
        if re.search(rf"\b{re.escape(loc)}\b", text):
            detected.add(loc)

    raw = {}
    for loc in detected:
        group = LOCATION_GROUP_MAP.get(loc, "unknown")
        raw[group] = raw.get(group, 0) + SUBCELLULAR_BASE_WEIGHTS[loc]

    if not raw:
        return "unknown", {}, {}

    total = sum(raw.values())
    scores = {k: round(v / total, 3) for k, v in raw.items()}
    dominant = max(scores, key=scores.get)

    return dominant, scores, {
        "ensembl_id": ensg,
        "detected_locations": sorted(detected),
        "multi_localizing": len(scores) > 1,
    }

# =========================
# DUB family inference
# =========================
def infer_dub_family(dub_gene_symbol):
    if dub_gene_symbol.startswith("USP"):
        return "USP"
    if dub_gene_symbol.startswith(("OTU", "OTUD")):
        return "OTU"
    if dub_gene_symbol.startswith("UCH"):
        return "UCH"
    if dub_gene_symbol.startswith("MINDY"):
        return "MINDY"
    if dub_gene_symbol.startswith("ATXN"):
        return "Josephin"
    return "Other"

# =========================
# UbiBrowser DUB fetch
# =========================
def fetch_ubibrowser_dubs(uniprot_ac):
    try:
        r = requests.get(
            f"http://ubibrowser.bio-it.cn/ubibrowser_v3/Home/Result/index/name/{uniprot_ac}/module/DUB/proteinType/noE3",
            timeout=30
        )
        r.raise_for_status()
        text = r.text
    except Exception:
        return {}

    def parse_js(var):
        m = re.search(rf"const\s+{var}\s*=\s*(\[[\s\S]*?\]);", text)
        if not m:
            return []
        try:
            return json.loads(re.sub(r",\s*([\}\]])", r"\1", m.group(1)))
        except Exception:
            return []

    dubs = {}
    for item in parse_js("knownTableData"):
        if item.get("targetGene"):
            dubs[item["targetGene"]] = {
                "confidence_score": 1.0,
                "level": "KNOWN"
            }

    for item in parse_js("predictedTableData"):
        g = item.get("targetGene")
        if g and g not in dubs:
            dubs[g] = {
                "confidence_score": float(item.get("SCORE", 0)),
                "level": "PREDICTED"
            }

    return dubs

# =========================
# FINAL INFERENCE (ENTRY POINT)
# =========================
def infer_local_dub_activity(substrate_uniprot_ac, e3_uniprot_ac=None):
    dubs = fetch_ubibrowser_dubs(substrate_uniprot_ac)
    if not dubs:
        return "low", {"source": "none"}

    substrate_gene = get_canonical_gene_symbol(substrate_uniprot_ac)
    sub_loc, _, _ = infer_subcellular_location_hpa_from_symbol(substrate_gene)

    e3_loc = None
    if e3_uniprot_ac:
        e3_gene = get_canonical_gene_symbol(e3_uniprot_ac)
        e3_loc, _, _ = infer_subcellular_location_hpa_from_symbol(e3_gene)

    enriched = {}
    for dub, info in dubs.items():
        family = infer_dub_family(dub)
        fam_score = DUB_FAMILY_BIAS.get(family, DUB_FAMILY_BIAS["Other"])["score"]

        dub_loc, _, _ = infer_subcellular_location_hpa_from_symbol(dub)

        match_score = 0.0
        if dub_loc == sub_loc:
            match_score += 0.5
        if e3_loc and dub_loc == e3_loc:
            match_score += 0.3
        if e3_loc and sub_loc == e3_loc:
            match_score += 0.2

        enriched[dub] = {
            **info,
            "family": family,
            "family_score": fam_score,
            "substrate_location": sub_loc,
            "dub_location": dub_loc,
            "e3_location": e3_loc,
            "colocalization_score": round(match_score, 3),
            "combined_dub_context_score": round(
                0.5 * fam_score
                + 0.3 * SUBCELLULAR_DUB_PRIORS.get(dub_loc, 0.3)
                + 0.2 * match_score,
                3
            ),
        }

    return "high", {"DUBs": enriched}
