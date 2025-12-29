import requests
import re

def is_multimeric_e3(uniprot_ac):
    url = (
        "https://reactome.org/ContentService/data/mapping/UniProt/"
        f"{uniprot_ac}/reactions?species=9606"
    )
    r = requests.get(url, headers={"accept": "application/json"}, timeout=20)
    r.raise_for_status()

    for entry in r.json():
        for name in entry.get("name", []):
            if ":" in name:
                return True
    return False
