# infinity_stone_3.py
import requests
import re

def is_multimeric_e3(uniprot_id):
   url = (
    f"https://reactome.org/ContentService/data/mapping/UniProt/"
    f"{UNIPROT_ID}/reactions?species={SPECIES}"
    )
    response = requests.get(url, headers={"accept": "application/json"})
    response.raise_for_status()
    data = response.json()

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

    is_multimeric = len(multimer_strings) > 0

    print("Reaction names containing multi-protein assemblies:\n")

    for m in sorted(multimer_strings):
        print("-", m)
    
    print("\nE3 multimeric:", is_multimeric)
    return False
