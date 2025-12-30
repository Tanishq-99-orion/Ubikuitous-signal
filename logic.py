###infinity stone 1

"""
E3 FLEXIBILITY SCORE (AUTONOMOUS, UbiBrowser-CORRECT)

Input:
- UniProt accession of E3 protein

Pipeline:
1. UniProt → gene symbol
2. Construct UbiBrowser term: GENE_HUMAN
3. UbiBrowser XML → enzyme_type (E3 class)
4. UniProt → domain coverage
5. MobiDB → disorder + pLDDT
6. Composite E3 flexibility score
"""

import requests
import xml.etree.ElementTree as ET

# ============================
# INPUT
# ============================
UNIPROT_AC = "Q96DX5" #E3 input

HEADERS_JSON = {"Accept": "application/json"}

# ============================
# URLS
# ============================
UNIPROT_URL = f"https://rest.uniprot.org/uniprotkb/{UNIPROT_AC}.json"
MOBIDB_URL  = f"https://mobidb.org/api/download_page?limit=5&acc={UNIPROT_AC}"

# ============================
# E3 CLASS → POLY PRIOR
# ============================
E3_POLYUB_SCORE = {
    "BTB_3": 0.85,
    "BTB_OTHER": 0.80,
    "CDC20": 0.90,
    "DWD": 0.80,
    "F-box": 0.85,
    "SOCS_VHL_BC-box": 0.80,
    "HECT": 0.45,
    "RING": 0.35,
    "UBOX": 0.30,
    "SINGLE_OTHER": 0.25
}

# ======================================================
# STEP 0: UniProt → construct UbiBrowser query term
# ======================================================
def get_ubibrowser_term_from_uniprot(uniprot_ac):
    r = requests.get(
        f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json",
        headers=HEADERS_JSON,
        timeout=20
    )
    r.raise_for_status()
    data = r.json()

    # Preferred: official gene symbol → GENE_HUMAN
    genes = data.get("genes", [])
    if genes:
        gene_symbol = genes[0].get("geneName", {}).get("value")
        if gene_symbol:
            return f"{gene_symbol}_HUMAN"

    # Fallback: UniProt entry name (already *_HUMAN)
    entry_name = data.get("uniProtkbId")
    if entry_name and entry_name.endswith("_HUMAN"):
        return entry_name

    return None


# ======================================================
# STEP 1: UbiBrowser → infer E3 class
# ======================================================
def get_e3_class_from_ubibrowser(ubibrowser_term):
    url = (
        "http://ubibrowser.ncpsb.org.cn/v2/home/API"
        f"?process=ubiquitination&type=enzyme&term={ubibrowser_term}"
    )

    r = requests.get(url, timeout=30)
    r.raise_for_status()

    root = ET.fromstring(r.text)

    e3_type_counts = {}

    for interaction in root.findall(".//Interaction"):
        species = interaction.findtext("species")
        etype = interaction.findtext("enzyme_type")

        if species == "H.sapiens" and etype:
            etype = etype.strip()
            e3_type_counts[etype] = e3_type_counts.get(etype, 0) + 1

    if not e3_type_counts:
        return None

    # Most frequent enzyme_type across interactions
    return max(e3_type_counts, key=e3_type_counts.get)


# ======================================================
# STEP 2: UniProt → domain coverage
# ======================================================
r = requests.get(UNIPROT_URL, headers=HEADERS_JSON, timeout=20)
r.raise_for_status()
up = r.json()

protein_length = up["sequence"]["length"]

KEEP_FEATURE_TYPES = {
    "Domain",
    "Repeat",
    "Region",
    "Zinc finger",
    "Coiled coil"
}

spans = []
features = []

for feat in up.get("features", []):
    if feat.get("type") not in KEEP_FEATURE_TYPES:
        continue

    loc = feat.get("location", {})
    start = loc.get("start", {}).get("value")
    end   = loc.get("end", {}).get("value")

    if start and end:
        spans.append((start, end))

    if feat.get("description"):
        features.append(f"{feat['type']}:{feat['description']}")

covered = set()
for s, e in spans:
    covered.update(range(s, e + 1))

domain_coverage = round(len(covered) / protein_length, 4)

# ======================================================
# STEP 3: Architecture score
# ======================================================
architecture_score = round(
    (1 - domain_coverage) * min(len(set(features)), 5) / 5,
    4
)

# ======================================================
# STEP 4: MobiDB → disorder + pLDDT
# ======================================================
r = requests.get(MOBIDB_URL, headers=HEADERS_JSON, timeout=30)
r.raise_for_status()
mobidb = r.json()["data"][0]

# ---- Disorder fraction
disorder_fraction = None
for key in [
    "prediction-disorder-priority",
    "homology-disorder-merge",
    "prediction-disorder-alphafold"
]:
    if key in mobidb and "content_fraction" in mobidb[key]:
        disorder_fraction = mobidb[key]["content_fraction"]
        break

if disorder_fraction is None:
    disorder_fraction = 0.0

# ---- pLDDT-based flexibility
plddt_block = mobidb.get("prediction-plddt-alphafold", {})
scores = plddt_block.get("scores", [])

frac_low = 0.0
frac_medium = 0.0

if scores:
    total = len(scores)
    low = sum(1 for s in scores if s < 0.5)
    medium = sum(1 for s in scores if 0.5 <= s < 0.7)

    frac_low = low / total
    frac_medium = medium / total

plddt_flex = round(
    1.0 * frac_low + 0.5 * frac_medium,
    4
)

# ======================================================
# STEP 5: E3 class inference
# ======================================================
ubibrowser_term = get_ubibrowser_term_from_uniprot(UNIPROT_AC)
if ubibrowser_term is None:
    raise RuntimeError("Failed to construct UbiBrowser term")

E3_CLASS = get_e3_class_from_ubibrowser(ubibrowser_term)
if E3_CLASS is None:
    E3_CLASS = "SINGLE_OTHER"

poly_prior = E3_POLYUB_SCORE.get(E3_CLASS, 0.3)

# ======================================================
# STEP 6: FINAL E3 FLEX SCORE
# ======================================================
E3_FLEX_SCORE = round(
    0.40 * poly_prior +
    0.25 * architecture_score +
    0.20 * disorder_fraction +
    0.15 * plddt_flex,
    4
)

# ======================================================
# OUTPUT
# ======================================================
print("\nE3 FLEXIBILITY SCORE (AUTONOMOUS)")
print("================================")
print("UniProt AC:", UNIPROT_AC)
print("UbiBrowser term:", ubibrowser_term)
print("E3 class (UbiBrowser):", E3_CLASS)
print("Polyubiquitination prior:", poly_prior)
print("Protein length:", protein_length)
print("Domain coverage:", domain_coverage)
print("Architecture score:", architecture_score)
print("Disorder fraction:", round(disorder_fraction, 4))
print("pLDDT frac_low:", round(frac_low, 4))
print("pLDDT frac_medium:", round(frac_medium, 4))
print("pLDDT flexibility:", plddt_flex)
print("FINAL E3_FLEX_SCORE:", E3_FLEX_SCORE)

###Infinity stone 2

import requests
import re

# =========================
# USER INPUT
# =========================
INPUT_PROTEIN = "Q96DX5"     # E3 protein (gene symbol)
SPECIES = 9606

# =========================
# API CONFIG
# =========================
REACTOME_URL = (
    f"https://reactome.org/ContentService/data/mapping/UniProt/"
    f"{INPUT_PROTEIN}/reactions?species={SPECIES}"
)

BIOGRID_URL = "https://webservice.thebiogrid.org/interactions/"
BIOGRID_ACCESS_KEY = "1ab875d99061e9ae48d022f346846e3b"

# =========================
# STEP 0: Gene symbol → UniProt accession
# =========================
def gene_symbol_to_uniprot(symbol, species=9606):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"(gene_exact:{symbol}) AND (organism_id:{species})",
        "fields": "accession",
        "format": "json",
        "size": 1
    }

    r = requests.get(url, params=params, timeout=15)
    if r.status_code != 200:
        return None

    results = r.json().get("results", [])
    if not results:
        return None

    return results[0].get("primaryAccession")


# =========================
# STEP 1: Reactome — extract E2s from reaction names
# =========================
def get_e2s_from_reactome():
    try:
        r = requests.get(REACTOME_URL, headers={"accept": "application/json"})
        r.raise_for_status()
        data = r.json()
    except Exception:
        return set()

    e2s = set()
    for entry in data:
        for name in entry.get("name", []):
            matches = re.findall(r"\bUBE2[A-Z0-9]+\b", name)
            e2s.update(matches)

    return e2s


# =========================
# STEP 2: BioGRID — fallback E2 extraction
# =========================
def get_e2s_from_biogrid(protein):
    params = {
        "geneList": protein,
        "format": "json",
        "accessKey": BIOGRID_ACCESS_KEY,
        "organism": SPECIES,
        "includeInteractors": "true",
        "interSpeciesExcluded": "true"
    }

    try:
        r = requests.get(BIOGRID_URL, params=params, timeout=30)
        r.raise_for_status()
        data = r.json()
    except Exception:
        return set()

    e2s = set()
    if isinstance(data, dict):
        for interaction in data.values():
            a = interaction.get("OFFICIAL_SYMBOL_A", "")
            b = interaction.get("OFFICIAL_SYMBOL_B", "")

            if a.upper() == protein.upper() and re.match(r"UBE2[A-Z0-9]+", b):
                e2s.add(b)
            if b.upper() == protein.upper() and re.match(r"UBE2[A-Z0-9]+", a):
                e2s.add(a)

    return e2s


# =========================
# STEP 3: UniProt — extract GO terms
# =========================
def get_go_terms_from_uniprot(uniprot_acc):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.json"
    r = requests.get(url, timeout=15)

    if r.status_code != 200:
        return []

    data = r.json()
    go_terms = []

    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GO":
            go_id = ref.get("id")
            desc = ""
            for prop in ref.get("properties", []):
                if prop.get("key") == "GoTerm":
                    desc = prop.get("value")
                    break
            if go_id and desc:
                go_terms.append((go_id, desc))

    return go_terms


# =========================
# STEP 4: E2 TYPE SCORING
# =========================
def score_e2_type(E2_gene, go_terms):
    score = {
        "priming": 0.0,
        "elongating": 0.0,
        "N_terminal": 0.0
    }

    # Canonical E2 priors
    canonical = {
        "UBE2D": "priming",
        "UBE2K": "elongating",
        "UBE2N": "elongating",
        "UBE2S": "elongating",
        "UBE2W": "N_terminal"
    }

    for prefix, etype in canonical.items():
        if E2_gene.startswith(prefix):
            score[etype] += 4.0

    for go_id, desc in go_terms:
        d = desc.lower()

        if "polyubiquitination" in d or "ubiquitin-ubiquitin" in d:
            score["elongating"] += 3.0

        elif any(k in d for k in [
            "k48-linked", "k63-linked", "k11-linked",
            "k6-linked", "k27-linked", "k29-linked", "k33-linked"
        ]):
            score["elongating"] += 3.0

        elif "linear polyubiquitination" in d or "n-end rule" in d:
            score["N_terminal"] += 3.5

        elif "monoubiquitination" in d:
            score["priming"] += 2.0

        elif "protein ubiquitination" in d:
            score["priming"] += 1.5

        elif "ubiquitin conjugating enzyme activity" in d:
            score["priming"] += 1.0

        elif "regulation of ubiquitin" in d:
            score["priming"] += 0.25

        elif any(x in d for x in [
            "deubiquitination", "binding", "reader",
            "desumoylation", "neddylation", "ufmylation"
        ]):
            continue

    sorted_scores = sorted(score.items(), key=lambda x: x[1], reverse=True)
    top, second = sorted_scores[0], sorted_scores[1]

        # --- revised dominance logic ---
    TOP_THRESHOLD = 3.0        # minimal evidence
    RATIO_THRESHOLD = 0.85     # how close scores must be to call "mixed"

    if top[1] < TOP_THRESHOLD:
        return "mixed", score

    if second[1] / top[1] >= RATIO_THRESHOLD:
        return "mixed", score

    return top[0], score

# =========================
# MAIN PIPELINE
# =========================
def main():
    print(f"\n🔍 Processing E3 protein: {INPUT_PROTEIN}")

    e2s = get_e2s_from_reactome()
    source = "Reactome"

    if not e2s:
        e2s = get_e2s_from_biogrid(INPUT_PROTEIN)
        source = "BioGRID"

    if not e2s:
        print("❌ No E2 enzymes found")
        return

    print(f"\n✅ E2 enzymes identified via {source}:\n")

    for e2 in sorted(e2s):
        print(f"🔹 E2: {e2}")

        uniprot_acc = gene_symbol_to_uniprot(e2)
        if not uniprot_acc:
            print("   ⚠ No UniProt accession found\n")
            continue

        go_terms = get_go_terms_from_uniprot(uniprot_acc)
        e2_type, score = score_e2_type(e2, go_terms)

        print(f"   🧬 UniProt: {uniprot_acc}")
        print(f"   🎯 E2 TYPE: {e2_type}")
        print(f"   📊 Score: {score}\n")


# =========================
# RUN
# =========================
if __name__ == "__main__":
    main()

###infinity stone 3

import requests
import re

# =========================
# USER INPUT
# =========================
UNIPROT_ID = "Q96DX5"   
SPECIES = 9606

# =========================
# REACTOME QUERY
# =========================
url = (
    f"https://reactome.org/ContentService/data/mapping/UniProt/"
    f"{UNIPROT_ID}/reactions?species={SPECIES}"
)

response = requests.get(url, headers={"accept": "application/json"})
response.raise_for_status()
data = response.json()

# =========================
# STEP 1: Extract multimer-like strings
# =========================
multimer_strings = set()

# Pattern:
# - contains ':' → interaction / complex
# - excludes gene names without modifiers
pattern = re.compile(r"\S+:\S+")

for entry in data:
    for name in entry.get("name", []):
        matches = pattern.findall(name)
        for m in matches:
            multimer_strings.add(m)

# =========================
# STEP 2: Decide multimeric or not
# =========================
is_multimeric = len(multimer_strings) > 0

# =========================
# OUTPUT
# =========================
print("Reaction names containing multi-protein assemblies:\n")

for m in sorted(multimer_strings):
    print("-", m)

print("\nE3 multimeric:", is_multimeric)

###infinity stone 4

import requests

# =========================
# USER INPUT
# =========================
SUBSTRATE_UNIPROT = "P04637"
SPECIES = 9606


# =========================
# Experimental ubiquitination sites
# (EBI Proteomics PTM API)
# =========================
def get_experimental_ub_sites_ebi(uniprot_acc, min_confidence={"Gold", "Silver"}):
    """
    Returns:
        ub_sites: set of ubiquitinated lysine positions (protein coordinates)
    """
    url = (
        "https://www.ebi.ac.uk/proteins/api/proteomics/ptm"
        f"?accession={uniprot_acc}&size=500"
    )

    headers = {"Accept": "application/json"}
    r = requests.get(url, headers=headers, timeout=30)

    if r.status_code != 200:
        return set()

    records = r.json()
    ub_sites = set()

    for record in records:
        for feature in record.get("features", []):
            if feature.get("type") != "PROTEOMICS_PTM":
                continue

            peptide_start = int(feature.get("begin", 0))

            for ptm in feature.get("ptms", []):
                if ptm.get("name") != "Ubiquitinylation":
                    continue

                ptm_pos = ptm.get("position")
                if not ptm_pos:
                    continue

                # Optional confidence filtering
                confidences = {
                    ref.get("properties", {}).get("Confidence score")
                    for ref in ptm.get("dbReferences", [])
                }

                if min_confidence and not (confidences & min_confidence):
                    continue

                protein_pos = peptide_start + int(ptm_pos) - 1
                ub_sites.add(protein_pos)

    return ub_sites


# =========================
# Accessible lysines estimator
# =========================
def get_accessible_lysines(uniprot_acc):
    """
    Estimate accessible lysines on a substrate protein.

    Returns:
        category: one | few | many
        details: dict
    """

    # -------------------------
    # 1. Experimental ubiquitination (MS-backed)
    # -------------------------
    ub_sites = get_experimental_ub_sites_ebi(uniprot_acc)
    num_ub_sites = len(ub_sites)

    # -------------------------
    # 2. UniProt sequence + disorder
    # -------------------------
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.json"
    r = requests.get(url, timeout=15)

    if r.status_code != 200:
        return "unknown", {}

    data = r.json()

    # Disordered regions
    disordered_ranges = []
    for feature in data.get("features", []):
        if feature.get("type") == "Region":
            if "disorder" in feature.get("description", "").lower():
                start = feature.get("location", {}).get("start", {}).get("value")
                end = feature.get("location", {}).get("end", {}).get("value")
                if start and end:
                    disordered_ranges.append((int(start), int(end)))

    # Lysine positions
    seq = data.get("sequence", {}).get("value", "")
    lys_positions = [i + 1 for i, aa in enumerate(seq) if aa == "K"]

    lys_in_idr = 0
    for k in lys_positions:
        for start, end in disordered_ranges:
            if start <= k <= end:
                lys_in_idr += 1
                break

    # -------------------------
    # 3. Encoding logic (your rule)
    # -------------------------
    if num_ub_sites >= 10:
        category = "many"
    elif num_ub_sites >= 3:
        category = "few"
    else:
        category = "one"

    details = {
        "experimentally_ubiquitinated_lysines": num_ub_sites,
        "lysines_in_disordered_regions": lys_in_idr,
        "total_lysines": len(lys_positions),
        "ubiquitinated_positions": sorted(ub_sites)
    }

    return category, details


# =========================
# MAIN
# =========================
def main():
    print(f"\n🔍 Substrate: {SUBSTRATE_UNIPROT}")

    category, details = get_accessible_lysines(SUBSTRATE_UNIPROT)

    print(f"\n📌 Accessible lysines category: **{category.upper()}**")
    print("📊 Details:")
    for k, v in details.items():
        print(f"   - {k}: {v}")


if __name__ == "__main__":
    main()

###infinity stone 5

import requests
import re
import json
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

GO_CONTEXT_MAP = {

    # =========================
    # ER-ASSOCIATED DEGRADATION
    # =========================
    "erad": {
        # Core ERAD identifiers
        "GO:0036503",  # ER-associated ubiquitin-dependent protein catabolic process
        "GO:0097466",  # ERAD pathway
        "GO:1904263",  # ERAD-mediated protein catabolic process

        # Supporting ubiquitin/proteasome logic (ER-linked)
        "GO:0045047",  # protein targeting to ER for degradation
        "GO:0016567",  # protein ubiquitination
        "GO:0043161",  # proteasome-mediated ubiquitin-dependent catabolic process

        # ER localization
        "GO:0005783",  # endoplasmic reticulum
        "GO:0031965",  # ER membrane
    },

    # =========================
    # EGAD (ENDOSOME / GOLGI)
    # =========================
    "egad": {
        # Endosomal / Golgi compartments (MANDATORY for EGAD)
        "GO:0005769",  # early endosome
        "GO:0005770",  # late endosome
        "GO:0032585",  # multivesicular body
        "GO:0005794",  # Golgi apparatus

        # Degradation machinery
        "GO:0006511",  # ubiquitin-dependent protein catabolic process
        "GO:0043161",  # proteasome-mediated ubiquitin-dependent catabolic process
        "GO:0016567",  # protein ubiquitination

        # Sorting / trafficking
        "GO:0006897",  # endocytosis
        "GO:0006898",  # receptor-mediated endocytosis
        "GO:0015031",  # protein transport
    },

    # =========================
    # SIGNALING (UBIQUITIN-RELEVANT)
    # =========================
    "signaling": {
        # General signaling
        "GO:0007165",  # signal transduction
        "GO:0023052",  # signaling

        # Immune / stress signaling
        "GO:0006955",  # immune response
        "GO:0007249",  # I-kappaB kinase/NF-kappaB signaling
        "GO:0043123",  # positive regulation of NF-kB signaling

        # MAPK / growth signaling
        "GO:0000165",  # MAPK cascade
        "GO:0070371",  # ERK1/2 cascade

        # Receptor signaling (often ubiquitin-controlled)
        "GO:0007169",  # RTK signaling
        "GO:0007186",  # GPCR signaling

        # Post-translational control nodes
        "GO:0018105",  # peptidyl-serine phosphorylation
        "GO:0006468",  # protein phosphorylation
    }
}

# ==================================================
# NEGATIVE EXCLUSION RULES (CRITICAL)
# ==================================================

GO_NEGATIVE_EXCLUSION = {

    # -------------------------
    # ERAD exclusions
    # -------------------------
    "erad": {
        "GO:0005769",  # early endosome
        "GO:0005770",  # late endosome
        "GO:0032585",  # multivesicular body
        "GO:0005794",  # Golgi apparatus
        "GO:0005886",  # plasma membrane
        "GO:0009897",  # external side of plasma membrane
    },

    # -------------------------
    # EGAD exclusions
    # -------------------------
    "egad": {
        "GO:0036503",  # ERAD pathway
        "GO:0097466",  # ER-associated degradation
        "GO:0045047",  # ER targeting for degradation
        "GO:0005783",  # ER
        "GO:0031965",  # ER membrane
    },

    # -------------------------
    # SIGNALING exclusions
    # -------------------------
    "signaling": {
        "GO:0006511",  # ubiquitin-dependent catabolism (alone ≠ signaling)
        "GO:0043161",  # proteasome-mediated degradation
        "GO:0005198",  # structural molecule activity
        "GO:0005737",  # cytosol (weak context alone)
    }
}

# ==================================================
# BIOGRID: INTERACTION COUNT
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
        data = r.json()
    except Exception:
        return 0

    return len(data)


# ==================================================
# REACTOME: PATHWAY → CURATED GO
# ==================================================

def reactome_pathway_to_go_accession(pathway_id):
    url = f"{REACTOME_BASE}/data/query/{pathway_id}"
    r = requests.get(url, headers={"accept": "application/json"}, timeout=20)
    r.raise_for_status()
    data = r.json()

    go_bp = data.get("goBiologicalProcess")
    if not go_bp:
        return None

    acc = go_bp.get("accession")
    if not acc:
        return None

    return f"GO:{acc}"

# ==================================================
# GO CONTEXT RESOLUTION (WITH EXCLUSIONS)
# ==================================================

CONTEXT_PRIORITY = ["erad", "egad", "signaling"]

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
            else:
                return ctx, "conflicted"

    return "other", "low"


def get_reactome_context_via_reactions(uniprot_id, debug=False):
    url = (
        f"{REACTOME_BASE}/data/mapping/UniProt/"
        f"{uniprot_id}/reactions?species=9606"
    )

    r = requests.get(url, headers={"accept": "application/json"}, timeout=20)
    r.raise_for_status()
    reactions = r.json()

    pathway_ids = set()

    for rxn in reactions:
        rxn_dbid = rxn.get("dbId")
        if not rxn_dbid:
            continue

        p_url = (
            f"{REACTOME_BASE}/search/pathways/of/"
            f"{rxn_dbid}"
            f"?includeInteractors=false&directlyInDiagram=false"
            f"&species={REACTOME_SPECIES.replace(' ', '%20')}"
        )

        r2 = requests.get(p_url, headers={"accept": "application/json"}, timeout=20)
        r2.raise_for_status()

        for p in r2.json():
            pid = p.get("stId") or p.get("dbId")
            if pid:
                pathway_ids.add(pid)

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


    return {
        "reaction_count": len(reactions),
        "pathway_count": len(pathway_ids),
        "go_terms": sorted(go_terms),
        "pathway_context": best_context
    }


# ==================================================
# DEGPRED PARSING
# ==================================================

def parse_degpred(uniprot_id):
    url = DEGPRED_URL.format(uniprot_id)
    r = requests.get(url, timeout=50)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")

    result = {
        "uniprot_id": uniprot_id,
        "sequence": None,
        "degpred_regions": [],
        "elm_regions": [],
        "disorder_regions": []
    }

    for script in soup.find_all("script"):
        if script.string and "FeatureViewer" in script.string:
            m = re.search(r'FeatureViewer\("([A-Z]+)"', script.string)
            if m:
                result["sequence"] = m.group(1)
            break

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
# RESIDENCE TIME INFERENCE (INTEGRATED)
# ==================================================

def overlaps(region, disorder_regions):
    return any(
        region["start"] <= d["end"] and region["end"] >= d["start"]
        for d in disorder_regions
    )


def analyze_residence_time(parsed, biogrid_count, reactome_ctx):
    degpred_regions = parsed["degpred_regions"]
    elm_regions = parsed["elm_regions"]
    disorder_regions = parsed["disorder_regions"]

    predicted_degrons = len(degpred_regions)
    degron_diversity = len(set(e.get("motif") for e in elm_regions))
    degrons_in_disorder = any(overlaps(d, disorder_regions) for d in degpred_regions)

    residence_score = 0

    if predicted_degrons >= 2:
        residence_score += 1

    if degron_diversity >= 2:
        residence_score += 1

    if degrons_in_disorder:
        residence_score += 1

    if biogrid_count >= 10:
        residence_score += 1

    if (
    reactome_ctx["pathway_context"] in {"egad", "erad"}
    and reactome_ctx.get("context_confidence") == "high"
):

        residence_score += 1

    residence_time = "long" if residence_score >= 3 else "short"

    return {
        "features": {
            "predicted_degrons": predicted_degrons,
            "degron_diversity": degron_diversity,
            "degrons_in_disorder": degrons_in_disorder,
            "biogrid_interactions": biogrid_count,
            "reactome_context": reactome_ctx["pathway_context"]
        },
        "residence_score": residence_score,
        "Substrate_residence_time": residence_time
    }


# ==================================================
# MASTER PIPELINE
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
        "biogrid": {"interaction_count": biogrid_count},
        "reactome": reactome_ctx,
        "residence_time_analysis": residence
    }


# =========================
# Example run
# =========================
if __name__ == "__main__":
    result = analyze_substrate("Q12888", "TP53")
    print(json.dumps(result, indent=2))

###infinity stone 6

# ==================================================
# GO HIERARCHY–BASED TRANSITION INFERENCE
# with caching + transition strength
# ==================================================

import requests
from functools import lru_cache

# -------------------------
# QuickGO API
# -------------------------
QUICKGO_ANCESTOR_URL = (
    "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
    "{go_id}/ancestors"
    "?relations=is_a,part_of,occurs_in,regulates"
)

HEADERS = {"Accept": "application/json"}

# -------------------------
# High-level GO anchors
# -------------------------

GO_ANCHORS = {

    # Stable functional state
    "signaling": {
        "GO:0007165",  # signal transduction
        "GO:0023052",  # signaling
        "GO:0006955",  # immune response
    },

    # Routing / transition processes
    "trafficking": {
        "GO:0006897",  # endocytosis
        "GO:0016192",  # vesicle-mediated transport
        "GO:0051179",  # localization
        "GO:0051234",  # establishment of localization
    },

    # EGAD commitment
    "egad_commitment": {
        "GO:0032585",  # multivesicular body
        "GO:0005769",  # early endosome
        "GO:0005770",  # late endosome
    },

    # ERAD commitment
    "erad_commitment": {
        "GO:0036503",  # ERAD
        "GO:0097466",  # ERAD pathway
        "GO:0045047",  # ER targeting for degradation
    },

    # Terminal degradation
    "terminal_degradation": {
        "GO:0006511",  # ubiquitin-dependent catabolism
        "GO:0043161",  # proteasome-mediated degradation
    }
}

# -------------------------
# Cached ancestor fetch
# -------------------------

@lru_cache(maxsize=2048)
def fetch_go_ancestors(go_id):
    """
    Fetch GO ancestors from QuickGO.
    Cached to avoid repeated API calls.
    """
    url = QUICKGO_ANCESTOR_URL.format(go_id=go_id.replace(":", "%3A"))

    r = requests.get(url, headers=HEADERS, timeout=20)
    r.raise_for_status()
    data = r.json()

    if not data.get("results"):
        return frozenset()

    ancestors = set(data["results"][0].get("ancestors", []))
    ancestors.add(go_id)

    return frozenset(ancestors)  # hashable for cache


# -------------------------
# Transition inference
# -------------------------

def infer_transition_from_go_hierarchy(go_id):
    """
    Infer transition context and strength from a SINGLE GO term
    using ontology hierarchy.
    """

    ancestors = fetch_go_ancestors(go_id)

    raw_hits = {
        category: ancestors & anchor_set
        for category, anchor_set in GO_ANCHORS.items()
    }

    # Serialize-friendly evidence
    hits = {k: sorted(v) for k, v in raw_hits.items() if v}

    # -------------------------
    # Transition strength score
    # -------------------------
    # Heuristic:
    # +1 per anchor category hit
    # +1 extra if commitment anchor present
    # +1 extra if terminal degradation present
    strength = 0

    for cat in hits:
        strength += 1

    if hits.get("egad_commitment") or hits.get("erad_commitment"):
        strength += 1

    if hits.get("terminal_degradation"):
        strength += 1

    # Cap for interpretability
    strength = min(strength, 5)

    # -------------------------
    # Transition classification
    # -------------------------

    if hits.get("signaling") and hits.get("trafficking"):
        transition = "signaling_to_egad"
        confidence = "high"

    elif hits.get("egad_commitment") and hits.get("terminal_degradation"):
        transition = "egad_to_degradation"
        confidence = "high"

    elif hits.get("erad_commitment"):
        transition = "signaling_or_qc_to_erad"
        confidence = "high"

    elif hits.get("trafficking"):
        transition = "routing_transition"
        confidence = "medium"

    else:
        transition = None
        confidence = "low"

    return {
        "go_term": go_id,
        "transition": transition,
        "transition_strength": strength,
        "confidence": confidence,
        "evidence": hits
    }

# ==================================================
# MAIN (STANDALONE TEST)
# ==================================================

if __name__ == "__main__":

    test_go_terms = [
        "GO:0006897",  # endocytosis
        "GO:0032585",  # multivesicular body
        "GO:0036503",  # ERAD
        "GO:0007165",  # signal transduction
        "GO:0006511",  # ubiquitin-dependent catabolism
    ]

    for go_id in test_go_terms:
        print("\n" + "=" * 60)
        print(f"GO TERM: {go_id}")

        try:
            result = infer_transition_from_go_hierarchy(go_id)
            print(json.dumps(result, indent=2))

        except Exception as e:
            print("Error:", e)

### infinity stone 7 
"""
SUBSTRATE-CENTRIC LOCAL DUB ACTIVITY INFERENCE
WITH DUB FAMILY BIAS + HPA SUBCELLULAR PRIORS
+ E3 / DUB SUBCELLULAR CO-LOCALIZATION MATCHING
(FINAL INTEGRATED VERSION — FIXED ENSEMBL RESOLUTION)
"""

import requests
import xml.etree.ElementTree as ET
import pandas as pd
import re
import json
from bs4 import BeautifulSoup

# =========================
# CONFIG
# =========================
BIOGRID_PATH = "/mnt/data/biogrid_E3_DUB_HUMAN_only.csv"

# =========================
# DUB FAMILY BIAS
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
# HPA SUBCELLULAR WEIGHTS
# =========================
SUBCELLULAR_BASE_WEIGHTS = {
    # nucleus
    "nucleoplasm": 1.0,
    "nuclear bodies": 0.9,
    "nuclear speckles": 0.85,
    "nucleoli": 0.8,
    "nuclear membrane": 0.75,
    "kinetochore": 0.7,
    "mitotic chromosome": 0.7,
    # endomembrane
    "endosomes": 1.0,
    "lysosomes": 0.9,
    "golgi apparatus": 0.7,
    "endoplasmic reticulum": 0.6,
    "vesicles": 0.6,
    "plasma membrane": 0.5,
    "lipid droplets": 0.4,
    "peroxisomes": 0.4,
    # cytosol
    "cytosol": 0.6,
    "cytoplasmic bodies": 0.6,
    "centrosome": 0.5,
    "mitotic spindle": 0.5,
    "midbody": 0.5,
    "microtubules": 0.45,
    "actin filaments": 0.4,
    "intermediate filaments": 0.35,
    "aggresome": 0.3,
    # mitochondria
    "mitochondria": 0.4,
}

LOCATION_GROUP_MAP = {
    # nucleus
    "nucleoplasm": "nucleus",
    "nuclear bodies": "nucleus",
    "nuclear speckles": "nucleus",
    "nucleoli": "nucleus",
    "nuclear membrane": "nucleus",
    "kinetochore": "nucleus",
    "mitotic chromosome": "nucleus",
    # endomembrane
    "endosomes": "endosome",
    "lysosomes": "endosome",
    "golgi apparatus": "golgi",
    "endoplasmic reticulum": "ERAD",
    "vesicles": "endosome",
    "plasma membrane": "endosome",
    "lipid droplets": "endomembrane",
    "peroxisomes": "endomembrane",
    # cytosol
    "cytosol": "cytosol",
    "cytoplasmic bodies": "cytosol",
    "centrosome": "cytosol",
    "mitotic spindle": "cytosol",
    "midbody": "cytosol",
    "microtubules": "cytosol",
    "actin filaments": "cytosol",
    "intermediate filaments": "cytosol",
    "aggresome": "cytosol",
    # mitochondria
    "mitochondria": "mitochondria",
}

SUBCELLULAR_DUB_PRIORS = {
    "nucleus": 0.8,
    "endosome": 0.85,
    "golgi": 0.6,
    "cytosol": 0.5,
    "proteasome": 0.3,
    "ERAD": 0.25,
    "mitochondria": 0.3,
    "endomembrane": 0.55,
    "unknown": 0.3,
}

# =========================
# UniProt helpers
# =========================
def get_canonical_gene_symbol(uniprot_ac):
    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json", timeout=15)
        r.raise_for_status()
        for gene in r.json().get("genes", []):
            if gene.get("geneName", {}).get("value"):
                return gene["geneName"]["value"]
    except requests.exceptions.RequestException:
        pass
    return None


def gene_symbol_to_uniprot(symbol, species=9606):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"(gene_exact:{symbol}) AND (organism_id:{species})",
        "fields": "accession",
        "format": "json",
        "size": 1,
    }
    try:
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        hits = r.json().get("results", [])
        if hits:
            return hits[0].get("primaryAccession")
    except requests.exceptions.RequestException:
        pass
    return None


# =========================
# FIXED: Gene symbol → Ensembl (no UniProt dependency)
# =========================
def gene_symbol_to_ensembl(symbol):
    """
    Uses Ensembl REST xrefs endpoint:
    /xrefs/symbol/homo_sapiens/{SYMBOL}
    """
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}"
    headers = {"Content-Type": "application/json"}

    try:
        r = requests.get(url, headers=headers, timeout=15)
        r.raise_for_status()
        data = r.json()
    except requests.exceptions.RequestException:
        return None

    for entry in data:
        if entry.get("type") == "gene" and entry.get("id", "").startswith("ENSG"):
            return entry["id"]

    return None


# =========================
# HPA subcellular inference
# =========================
def infer_subcellular_location_hpa_from_symbol(gene_symbol):
    ensg = gene_symbol_to_ensembl(gene_symbol)
    if not ensg:
        return "unknown", {}, {}

    try:
        r = requests.get(f"https://www.proteinatlas.org/{ensg}/subcellular", timeout=20)
        r.raise_for_status()
        soup = BeautifulSoup(r.text, "html.parser")
        text = soup.get_text(" ").lower()
    except requests.exceptions.RequestException:
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
# DUB family
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
# UbiBrowser v3 DUB fetch
# =========================
def fetch_ubibrowser_dubs(uniprot_ac):
    try:
        r = requests.get(
            f"http://ubibrowser.bio-it.cn/ubibrowser_v3/Home/Result/index/name/{uniprot_ac}/module/DUB/proteinType/noE3",
            timeout=30,
        )
        r.raise_for_status()
        text = r.text
    except requests.exceptions.RequestException:
        return {}

    def parse_js(var):
        m = re.search(rf"const\s+{var}\s*=\s*(\[[\s\S]*?\]);", text)
        if not m:
            return []
        js = re.sub(r",\s*([\}\]])", r"\1", m.group(1))
        try:
            return json.loads(js)
        except Exception:
            return []

    dubs = {}
    for item in parse_js("knownTableData"):
        if item.get("targetGene"):
            dubs[item["targetGene"]] = {"confidence_score": 1.0, "level": "KNOWN"}

    for item in parse_js("predictedTableData"):
        g = item.get("targetGene")
        if g and g not in dubs:
            dubs[g] = {"confidence_score": float(item.get("SCORE", 0)), "level": "PREDICTED"}

    return dubs


# =========================
# Enrich DUB + E3 context
# =========================
def enrich_dub_context(dubs_dict, substrate_uniprot_ac, e3_uniprot_ac=None):
    substrate_gene = get_canonical_gene_symbol(substrate_uniprot_ac)
    sub_loc, _, _ = infer_subcellular_location_hpa_from_symbol(substrate_gene)

    e3_loc = None
    if e3_uniprot_ac:
        e3_gene = get_canonical_gene_symbol(e3_uniprot_ac)
        e3_loc, _, _ = infer_subcellular_location_hpa_from_symbol(e3_gene)

    enriched = {}
    for dub_gene_symbol, info in dubs_dict.items():
        family = infer_dub_family(dub_gene_symbol)
        fam_score = DUB_FAMILY_BIAS.get(family, DUB_FAMILY_BIAS["Other"])["score"]

        dub_loc, _, _ = infer_subcellular_location_hpa_from_symbol(dub_gene_symbol)

        match_score = 0.0
        if dub_loc == sub_loc:
            match_score += 0.5
        if e3_loc and dub_loc == e3_loc:
            match_score += 0.3
        if e3_loc and sub_loc == e3_loc:
            match_score += 0.2

        enriched[dub_gene_symbol] = {
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
                3,
            ),
        }

    return enriched


# =========================
# FINAL INFERENCE
# =========================
def infer_local_dub_activity(substrate_uniprot_ac, e3_uniprot_ac=None):
    dubs = fetch_ubibrowser_dubs(substrate_uniprot_ac)
    if not dubs:
        return "low", {"source": "none"}

    enriched = enrich_dub_context(dubs, substrate_uniprot_ac, e3_uniprot_ac)
    if not enriched:
        return "low", {"source": "no_hpa_mapping"}

    return "high", {"DUBs": enriched}


# =========================
# MAIN
# =========================
if __name__ == "__main__":

    substrate_uniprot_ac = "P04637"  # TP53
    e3_uniprot_ac = "Q00987"         # MDM2

    state, evidence = infer_local_dub_activity(
        substrate_uniprot_ac=substrate_uniprot_ac,
        e3_uniprot_ac=e3_uniprot_ac
    )

    print("\n==============================")
    print("LOCAL DUB ACTIVITY INFERENCE")
    print("==============================")
    print("Substrate UniProt:", substrate_uniprot_ac)
    print("E3 UniProt:", e3_uniprot_ac)
    print("Local_DUB_activity:", state)

    if state == "high":
        print("\n--- DUB CONTEXT DETAILS ---")
        for dub, info in evidence["DUBs"].items():
            print(f"\nDUB: {dub}")
            for k, v in info.items():
                print(f"  {k:25s}: {v}")

    print("==============================\n")

###infinity stone 8

"""
BIOCHEMICAL LINKAGE GRAMMAR LAYER
Rule-based override for non-canonical ubiquitin chain types
(INPUTS: substrate UniProt AC, E3 ligase gene symbol)
Designed to sit ON TOP of your systems-level model
"""

# =========================
# LINKAGE RULE DEFINITIONS
# =========================

# E2-driven linkage specificity (hard biochemical bias)
E2_LINKAGE_RULES = {
    "UBE2S":  {"K11": 0.95},
    "UBE2N":  {"K63": 0.9},
    "UBE2V1": {"K63": 0.9},
    "UBE2V2": {"K63": 0.9},
    "UBE2D":  {"K48": 0.5, "K63": 0.5},  # promiscuous
}

# E3-driven linkage enforcement (instructional E3s)
E3_LINKAGE_RULES = {
    "APC/C":      {"K11": 1.0},
    "ANAPC":      {"K11": 1.0},
    "BRCA1":      {"K6": 1.0},
    "BARD1":      {"K6": 1.0},
    "HOIP":       {"M1": 1.0},
    "RNF31":      {"M1": 1.0},
    "TRIM13":     {"K11": 0.8},
    "TRIM25":     {"K63": 0.7, "K48": 0.3},
}

# DUBs that actively erase specific chains (negative pressure)
DUB_ANTI_LINKAGE = {
    "OTULIN": {"M1": -1.0},
    "Cezanne": {"K11": -0.8},
    "TRABID": {"K29": -0.8, "K33": -0.8},
    "CYLD": {"K63": -0.6},
}

# =========================
# HELPERS
# =========================
def normalize_gene(g):
    return g.upper().replace("-", "").replace("_", "")


def resolve_e3_family(e3_gene):
    """
    Coarse E3 family resolution from gene name
    """
    g = normalize_gene(e3_gene)
    if g.startswith("ANAPC"):
        return "APC/C"
    if g in {"BRCA1", "BARD1"}:
        return g
    if g in {"RNF31", "HOIP"}:
        return "HOIP"
    return g


# =========================
# MAIN BIOCHEMICAL LAYER
# =========================
def biochemical_linkage_override(
    substrate_ac,
    e3_gene,
    e2_gene=None,
    local_dubs=None,
):
    """
    Decide non-canonical ubiquitin linkage types using biochemical rules.

    Inputs:
        substrate_ac : UniProt accession (kept for extensibility)
        e3_gene      : E3 ligase gene symbol
        e2_gene      : optional E2 gene symbol
        local_dubs   : dict from your DUB module (optional)

    Returns:
        dict with enforced / biased linkage outcomes
    """

    linkage_scores = {}

    # -------------------------
    # 1. E3 instructional rules
    # -------------------------
    e3_key = resolve_e3_family(e3_gene)
    if e3_key in E3_LINKAGE_RULES:
        for link, score in E3_LINKAGE_RULES[e3_key].items():
            linkage_scores[link] = linkage_scores.get(link, 0) + score

    # -------------------------
    # 2. E2 catalytic bias
    # -------------------------
    if e2_gene:
        e2_key = normalize_gene(e2_gene)
        if e2_key in E2_LINKAGE_RULES:
            for link, score in E2_LINKAGE_RULES[e2_key].items():
                linkage_scores[link] = linkage_scores.get(link, 0) + score

    # -------------------------
    # 3. DUB anti-linkage pressure
    # -------------------------
    if local_dubs:
        for dub in local_dubs:
            dub_key = normalize_gene(dub)
            if dub_key in DUB_ANTI_LINKAGE:
                for link, penalty in DUB_ANTI_LINKAGE[dub_key].items():
                    linkage_scores[link] = linkage_scores.get(link, 0) + penalty

    # -------------------------
    # 4. Normalize & classify
    # -------------------------
    if not linkage_scores:
        return {
            "mode": "system_level_only",
            "allowed_linkages": ["K48", "K63", "mixed"],
            "scores": {},
        }

    # Clamp scores
    linkage_scores = {
        k: round(max(min(v, 1.0), -1.0), 3)
        for k, v in linkage_scores.items()
    }

    dominant = max(linkage_scores, key=linkage_scores.get)

    return {
        "mode": "biochemical_override",
        "dominant_linkage": dominant,
        "scores": linkage_scores,
    }


# =========================
# EXAMPLE
# =========================
if __name__ == "__main__":
    result = biochemical_linkage_override(
        substrate_ac="P04637",
        e3_gene="ANAPC1",
        e2_gene="UBE2S",
        local_dubs={"USP7": {}, "OTUB1": {}},
    )

    print("Biochemical linkage decision:")
    print(result)

### infinity stone 9

"""
TEMPORAL / CELL-STATE GATING MODULE
(EXPLICIT CELL STATE PROVIDED BY USER / EXPERIMENT)
Hard permission layer for ubiquitin linkage prediction
"""

# =========================
# ALLOWED CELL STATES
# =========================
CELL_STATES = {
    "quiescent",
    "cycling",
    "DNA_damage",
    "immune_active",
}

# =========================
# TEMPORAL GATING RULES
# =========================
TEMPORAL_GATES = {
    "cycling": {
        "allow": {"K11", "K48"},
        "block": {"K6", "K27"},
    },
    "DNA_damage": {
        "allow": {"K6", "K63", "mixed"},
        "block": {"K11"},
    },
    "immune_active": {
        "allow": {"K27", "K63"},
        "block": {"K11", "K6"},
    },
    "quiescent": {
        "allow": {"K48", "K63"},
        "block": {"K11", "K6", "K27"},
    },
}

# =========================
# USER CELL-STATE INPUT
# =========================
def get_user_cell_state():
    """
    Force user to explicitly choose a valid cell state.
    """
    print("\nSelect experimental / cellular state:")
    for state in sorted(CELL_STATES):
        print(f" - {state}")

    while True:
        choice = input("\nEnter cell state: ").strip()

        if choice in CELL_STATES:
            return choice

        print(
            f"❌ Invalid cell state '{choice}'. "
            f"Choose one of: {sorted(CELL_STATES)}"
        )


# =========================
# APPLY TEMPORAL GATE
# =========================
def apply_temporal_gate_explicit(linkage_scores, cell_state):
    """
    Apply temporal gating using an explicitly provided cell state.
    """

    gate = TEMPORAL_GATES[cell_state]

    allowed = gate["allow"]
    blocked = gate["block"]

    # Remove forbidden linkages
    gated = {
        k: v for k, v in linkage_scores.items()
        if k in allowed
    }

    if not gated:
        return {}, {
            "cell_state": cell_state,
            "gate_applied": True,
            "allowed": sorted(allowed),
            "blocked": sorted(blocked),
            "note": "All candidate linkages blocked by temporal gate",
        }

    # Renormalize scores
    total = sum(gated.values())
    if total > 0:
        gated = {k: round(v / total, 3) for k, v in gated.items()}

    return gated, {
        "cell_state": cell_state,
        "gate_applied": True,
        "allowed": sorted(allowed),
        "blocked": sorted(blocked),
    }


# =========================
# EXAMPLE RUN
# =========================
if __name__ == "__main__":

    raw_linkage_scores = {
        "K11": 0.5,
        "K48": 0.3,
        "K63": 0.2,
        "K6": 0.1,
    }

    # 🔒 HARD USER INPUT GATE
    cell_state = get_user_cell_state()

    gated_scores, info = apply_temporal_gate_explicit(
        raw_linkage_scores,
        cell_state
    )

    print("\n[INFO] Cell state:", cell_state)
    print("[INFO] Gate info:", info)
    print("[RESULT] Gated linkage scores:", gated_scores)

"""
INFINITY STONE 10
=================
Cell-type + ubiquitin-linkage → process bias mapper (BN-friendly)

Inputs:
- *_cell_type_priors.csv files in working directory
  (lysosome, proteasome, autophagy, escrt, recycling, deubiquitination)

Interactive notebook usage:
- User selects one of 15 cell types
- User selects ubiquitin linkage (including MONO)
- Outputs ranked posterior probabilities over processes
"""

import os
import glob
import argparse
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
    "Stem and proliferating cells"
]

UBIQUITIN_LINKAGES = [
    "MONO",
    "K48", "K63", "K11", "K29",
    "K6", "K27", "K33",
    "M1"
]

# =========================================================
# LINKAGE → PROCESS BIASES  P(linkage | process)
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

# Normalize likelihoods
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
    "M1":  0.2
}

EPSILON_PRIOR = 1e-6  # used if a cell type is missing for a process

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


def select_from_list(options, title):
    print(f"\nSelect {title}:")
    for i, opt in enumerate(options):
        print(f"  [{i}] {opt}")
    idx = int(input("Enter index: ").strip())
    return options[idx]

# =========================================================
# CORE MAPPING FUNCTION
# =========================================================

def map_linkage_to_process(
    cell_type: str,
    linkage: str,
    priors_dict: Dict[str, pd.Series],
    dub_process_name: str = "deubiquitination"
) -> Tuple[pd.DataFrame, Dict]:

    linkage = linkage.upper()
    processes = sorted(priors_dict.keys())

    # P(process | cell)
    p_cell = {}
    for p in processes:
        s = priors_dict[p]
        p_cell[p] = float(s.loc[cell_type]) if cell_type in s.index else EPSILON_PRIOR

    # P(linkage | process)
    p_l_given_proc = {
        p: LINKAGE_LIKELIHOOD.get(linkage, {}).get(p, 1e-6)
        for p in processes
    }

    # DUB prior
    dub_series = priors_dict.get(dub_process_name, pd.Series())
    dub_prior = float(dub_series.loc[cell_type]) if cell_type in dub_series.index else 0.0
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
            "Unnormalized_score": unnorm
        })

    df = pd.DataFrame(rows)

    if df.empty or "Unnormalized_score" not in df.columns:
           return (
             pd.DataFrame(
                columns=[
                "Process",
                "Posterior",
                "Prior_given_cell",
                "P_linkage_given_process",
                "Survival_factor"
                 ]
             ),
           {
            "cell_type": cell_type,
            "linkage": linkage,
            "error": "no_process_priors_loaded"
           }
           )

    total = df["Unnormalized_score"].sum()

    df["Posterior"] = (
        df["Unnormalized_score"] / total
        if total > 0
        else df["Prior_given_cell"] / df["Prior_given_cell"].sum()
    )

    df = df.sort_values("Posterior", ascending=False).reset_index(drop=True)

    metadata = {
        "cell_type": cell_type,
        "linkage": linkage,
        "dub_prior": dub_prior,
        "editability": editability
    }

    return df, metadata

# =========================================================
# MAIN (NOTEBOOK MODE)
# =========================================================

if __name__ == "__main__":

    prior_files = find_prior_files()
    if not prior_files:
        raise FileNotFoundError("No *_cell_type_priors.csv files found")

    priors_dict = load_priors(prior_files)

    CELL_TYPE = select_from_list(CELL_TYPE_CHOICES, "cell type")
    LINKAGE = select_from_list(UBIQUITIN_LINKAGES, "ubiquitin linkage")

    df, meta = map_linkage_to_process(CELL_TYPE, LINKAGE, priors_dict)

    print("\n=== INFERENCE RESULT ===")
    print(f"Cell type : {meta['cell_type']}")
    print(f"Linkage   : {meta['linkage']}")
    print(f"DUB prior : {meta['dub_prior']:.3f}\n")

    display(
        df[[
            "Process",
            "Posterior",
            "Prior_given_cell",
            "P_linkage_given_process",
            "Survival_factor"
        ]]
    )

