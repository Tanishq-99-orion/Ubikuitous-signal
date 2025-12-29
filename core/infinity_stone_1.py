import requests
import xml.etree.ElementTree as ET

HEADERS_JSON = {"Accept": "application/json"}

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

def get_ubibrowser_term_from_uniprot(uniprot_ac):
    r = requests.get(
        f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json",
        headers=HEADERS_JSON,
        timeout=20
    )
    r.raise_for_status()
    data = r.json()

    genes = data.get("genes", [])
    if genes:
        g = genes[0].get("geneName", {}).get("value")
        if g:
            return f"{g}_HUMAN"

    entry = data.get("uniProtkbId")
    if entry and entry.endswith("_HUMAN"):
        return entry

    return None

def get_e3_class_from_ubibrowser(term):
    url = (
        "http://ubibrowser.ncpsb.org.cn/v2/home/API"
        f"?process=ubiquitination&type=enzyme&term={term}"
    )
    r = requests.get(url, timeout=30)
    r.raise_for_status()

    root = ET.fromstring(r.text)
    counts = {}

    for it in root.findall(".//Interaction"):
        if it.findtext("species") == "H.sapiens":
            etype = it.findtext("enzyme_type")
            if etype:
                counts[etype] = counts.get(etype, 0) + 1

    return max(counts, key=counts.get) if counts else None
