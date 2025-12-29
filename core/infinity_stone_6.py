import requests
from functools import lru_cache

QUICKGO_ANCESTOR_URL = (
    "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
    "{go_id}/ancestors"
    "?relations=is_a,part_of,occurs_in,regulates"
)

HEADERS = {"Accept": "application/json"}

GO_ANCHORS = {
    "signaling": {
        "GO:0007165",
        "GO:0023052",
        "GO:0006955",
    },
    "trafficking": {
        "GO:0006897",
        "GO:0016192",
        "GO:0051179",
        "GO:0051234",
    },
    "egad_commitment": {
        "GO:0032585",
        "GO:0005769",
        "GO:0005770",
    },
    "erad_commitment": {
        "GO:0036503",
        "GO:0097466",
        "GO:0045047",
    },
    "terminal_degradation": {
        "GO:0006511",
        "GO:0043161",
    }
}

@lru_cache(maxsize=2048)
def fetch_go_ancestors(go_id):
    url = QUICKGO_ANCESTOR_URL.format(go_id=go_id.replace(":", "%3A"))
    r = requests.get(url, headers=HEADERS, timeout=20)
    r.raise_for_status()

    data = r.json()
    if not data.get("results"):
        return frozenset()

    ancestors = set(data["results"][0].get("ancestors", []))
    ancestors.add(go_id)
    return frozenset(ancestors)


def infer_transition_from_go_hierarchy(go_id):
    ancestors = fetch_go_ancestors(go_id)

    hits = {
        category: sorted(ancestors & anchor_set)
        for category, anchor_set in GO_ANCHORS.items()
        if ancestors & anchor_set
    }

    strength = len(hits)
    if hits.get("egad_commitment") or hits.get("erad_commitment"):
        strength += 1
    if hits.get("terminal_degradation"):
        strength += 1

    strength = min(strength, 5)

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
