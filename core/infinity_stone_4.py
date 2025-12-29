import requests

# =========================
# Experimental ubiquitination sites
# (EBI Proteomics PTM API)
# =========================
def get_experimental_ub_sites_ebi(uniprot_acc, min_confidence={"Gold", "Silver"}):
    """
    Returns:
        ub_sites: set of ubiquitinated lysine positions
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

    ub_sites = get_experimental_ub_sites_ebi(uniprot_acc)
    num_ub_sites = len(ub_sites)

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.json"
    r = requests.get(url, timeout=15)

    if r.status_code != 200:
        return "unknown", {}

    data = r.json()

    disordered_ranges = []
    for feature in data.get("features", []):
        if feature.get("type") == "Region":
            if "disorder" in feature.get("description", "").lower():
                start = feature.get("location", {}).get("start", {}).get("value")
                end = feature.get("location", {}).get("end", {}).get("value")
                if start and end:
                    disordered_ranges.append((int(start), int(end)))

    seq = data.get("sequence", {}).get("value", "")
    lys_positions = [i + 1 for i, aa in enumerate(seq) if aa == "K"]

    lys_in_idr = 0
    for k in lys_positions:
        for start, end in disordered_ranges:
            if start <= k <= end:
                lys_in_idr += 1
                break

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
