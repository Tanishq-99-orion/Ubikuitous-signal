import os
import re
import requests
from dotenv import load_dotenv

load_dotenv()

SPECIES = 9606
BIOGRID_URL = "https://webservice.thebiogrid.org/interactions"
BIOGRID_ACCESS_KEY = os.getenv("BIOGRID_ACCESS_KEY")

def gene_symbol_to_uniprot(symbol, species=SPECIES):
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
    hits = r.json().get("results", [])
    return hits[0]["primaryAccession"] if hits else None

def get_e2s_from_reactome(uniprot_ac=None):
    if not uniprot_ac:
        return set()
    url = (
        "https://reactome.org/ContentService/data/mapping/UniProt/"
        f"{uniprot_ac}/reactions?species={SPECIES}"
    )
    try:
        r = requests.get(url, headers={"accept": "application/json"}, timeout=20)
        r.raise_for_status()
    except Exception:
        return set()

    e2s = set()
    for entry in r.json():
        for name in entry.get("name", []):
            e2s.update(re.findall(r"\bUBE2[A-Z0-9]+\b", name))
    return e2s

def get_e2s_from_biogrid(gene):
    if not gene or not BIOGRID_ACCESS_KEY:
        return set()

    params = {
        "geneList": gene,
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
    for it in data.values():
        a = it.get("OFFICIAL_SYMBOL_A", "")
        b = it.get("OFFICIAL_SYMBOL_B", "")
        if a == gene and b.startswith("UBE2"):
            e2s.add(b)
        if b == gene and a.startswith("UBE2"):
            e2s.add(a)

    return e2s

def get_go_terms_from_uniprot(uniprot_ac):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json"
    r = requests.get(url, timeout=15)
    if r.status_code != 200:
        return []

    gos = []
    for ref in r.json().get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GO":
            for p in ref.get("properties", []):
                if p.get("key") == "GoTerm":
                    gos.append((ref["id"], p["value"]))
    return gos

def score_e2_type(e2_gene, go_terms):
    score = {"priming": 0, "elongating": 0, "N_terminal": 0}

    if e2_gene.startswith("UBE2S"):
        score["elongating"] += 4
    if e2_gene.startswith("UBE2N"):
        score["elongating"] += 4
    if e2_gene.startswith("UBE2W"):
        score["N_terminal"] += 4

    for _, desc in go_terms:
        d = desc.lower()
        if "polyubiquitination" in d:
            score["elongating"] += 2
        elif "monoubiquitination" in d:
            score["priming"] += 2

    top = max(score, key=score.get)
    return top, score
