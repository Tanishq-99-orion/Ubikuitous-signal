# infinity_stone_2.py
import os
import re
import requests
from dotenv import load_dotenv

load_dotenv()

BIOGRID_ACCESS_KEY = os.getenv("BIOGRID_ACCESS_KEY")
SPECIES = 9606

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
    res = r.json().get("results", [])
    return res[0]["primaryAccession"] if res else None

def get_e2s_from_reactome():
    # paste your Reactome E2 extraction logic here
    return set()

def get_e2s_from_biogrid(gene):
    # paste your BioGRID E2 extraction logic here
    return set()

def get_go_terms_from_uniprot(uniprot_ac):
    # paste your GO extraction logic here
    return []

def score_e2_type(e2_gene, go_terms):
    # paste your scoring logic here
    return "mixed", {}
