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
# APPLY TEMPORAL GATE
# =========================
def apply_temporal_gate_explicit(linkage_scores, cell_state):
    """
    Apply temporal gating using an explicitly provided cell state.

    Inputs:
        linkage_scores : dict[str, float]
        cell_state     : str (must be in CELL_STATES)

    Returns:
        gated_scores : dict
        gate_info    : dict (metadata)
    """

    if cell_state not in CELL_STATES:
        return linkage_scores, {
            "cell_state": cell_state,
            "gate_applied": False,
            "reason": "invalid_cell_state",
        }

    gate = TEMPORAL_GATES[cell_state]
    allowed = gate["allow"]
    blocked = gate["block"]

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
            "note": "all_candidate_linkages_blocked",
        }

    total = sum(gated.values())
    if total > 0:
        gated = {k: round(v / total, 3) for k, v in gated.items()}

    return gated, {
        "cell_state": cell_state,
        "gate_applied": True,
        "allowed": sorted(allowed),
        "blocked": sorted(blocked),
    }
