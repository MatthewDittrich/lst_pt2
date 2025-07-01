import awkward as ak
import numpy as np

from tools import pt2_methods


class LST_OBJECT:

    def __init__(self, data: dict):
        if not data:
            raise ValueError("Cannot create LST_OBJECT with empty data.")
        # Get reference attr
        ref_key = next(iter(data))
        ref_array = data[ref_key]
        # Check other attr
        for k, v in data.items():
            if len(v) != len(ref_array):
                raise ValueError(f"Outer length mismatch: '{k}' has length {len(v)}, expected {len(ref_array)}")
            if not ak.all(ak.num(v) == ak.num(ref_array)):
                raise ValueError(f"Inner length mismatch: '{k}' has inconsistent per-event structure compared to '{ref_key}'")
        self._data = data
    
    def __getattr__(self, key):
        try:
            return self._data[key]
        except KeyError:
            raise AttributeError(f"'LST_OBJECT' has no attribute '{key}'")

    def __getitem__(self, key):
        if isinstance(key, str):
            return self._data[key]
        elif isinstance(key, ak.Array):
            # Assume this is a mask; apply it to all fields
            masked_data = {k: v[key] for k, v in self._data.items()}
            return LST_OBJECT(masked_data)
        else:
            raise TypeError(f"Unsupported key type: {type(key)}")

    def __setitem__(self, key, value):
        if not self._data:
            raise ValueError("Cannot add a field to an empty LST_OBJECT.")
        # Reference attr to compare against
        ref_key = next(iter(self._data))
        ref_array = self._data[ref_key]
        # Outer length check
        if len(value) != len(ref_array):
            raise ValueError(f"Outer length mismatch: '{key}' has length {len(value)} but expected {len(ref_array)}")
        # Inner length check (per event)
        if not ak.all(ak.num(value) == ak.num(ref_array)):
            raise ValueError(f"Inner length mismatch: '{key}' has inconsistent per-event structure compared to '{ref_key}'")
        self._data[key] = value

    def __repr__(self):
        return f"LST_OBJECT(fields={list(self._data.keys())})"


def get_ls(events):
    
    # Define the object and the basics
    ls = LST_OBJECT({
        "pt": events["ls_pt"],
        "eta": events["ls_eta"],
        "phi": events["ls_phi"],
        "simIdx": events["ls_simIdx"],
    })
    ls["isfake"]  = (events.ls_isFake == 1)
    ls["is_from_pt5"] = pt2_methods.ls_masks(events, "pt5")
    ls["is_from_pt3"] = pt2_methods.ls_masks(events, "pt3")
    ls["is_from_t5"] = pt2_methods.ls_masks(events, "t5")
    ls["is_used"] = (ls.is_from_pt5 | ls.is_from_t5 | ls.is_from_pt3)

    return ls


def get_pls(events):

    # Define the object and the basics
    pls = LST_OBJECT({
        "pt": events["pls_pt"],
        "eta": events["pls_eta"],
        "phi": events["pls_phi"],
        "simIdx": events["pls_simIdx"],
    })
    pls["isfake"] = (events.pls_isFake == 1)
    pls["is_from_pt5"] = pt2_methods.pls_masks(events, "pt5")
    pls["is_from_tcpls"] = pt2_methods.pls_masks(events, "pLS")
    pls["is_from_pt3"] = pt2_methods.pls_masks(events, "pt3")
    pls["is_used"] = (pls.is_from_pt5 | pls.is_from_pt3)

    return pls


def get_ideal_pt2(pls, ls):
    filtered_pls = pls[pls.matches_unused_ls]
    filtered_ls = ls[ls.matches_unused_pls]

    filtered_pls_simIdx = filtered_pls.simIdx
    filtered_ls_simIdx = filtered_ls.simIdx

    pls_idx = []
    ls_idx = []

    # Loop per event (not ideal but works)
    for pls_sim, ls_sim in zip(filtered_pls_simIdx, filtered_ls_simIdx):
        # Build a dict of ls_simIdx vals and respective ids
        sim_to_ls = {}
        for j, sim in enumerate(ls_sim):
            sim_to_ls.setdefault(sim, []).append(j)

        event_pls_idx = []
        event_ls_idx = []

        # Loop over each pls, match to all ls entries via dict
        for i, sim in enumerate(pls_sim):
            if sim in sim_to_ls:
                matched_ls_indices = sim_to_ls[sim]
                event_pls_idx.extend([i] * len(matched_ls_indices))  # same pls index repeated
                event_ls_idx.extend(matched_ls_indices)

        pls_idx.append(event_pls_idx)
        ls_idx.append(event_ls_idx)

    pls_idx = ak.Array(pls_idx)
    ls_idx  = ak.Array(ls_idx)

    pt2 = LST_OBJECT({
        "pls_pt": filtered_pls.pt[pls_idx],
        "pls_eta": filtered_pls.eta[pls_idx],
        "pls_phi": filtered_pls.phi[pls_idx],
        "ls_pt": filtered_ls.pt[ls_idx],
        "ls_eta": filtered_ls.eta[ls_idx],
        "ls_phi": filtered_ls.phi[ls_idx],
        "simIdx": filtered_ls.simIdx[ls_idx],  # or pls.simIdx[pls_idx], same here
    })

    # Compute some values
    pt2["delta_pt"] = pt2.pls_pt - pt2.ls_pt
    pt2["delta_eta"] = pt2.pls_eta - pt2.ls_eta
    delta_phi_raw = pt2.pls_phi - pt2.ls_phi
    pt2["delta_phi"] = np.arctan2(np.sin(delta_phi_raw), np.cos(delta_phi_raw))
    pt2["delta_R"] = np.sqrt(pt2.delta_eta**2 + pt2.delta_phi**2)

    # Right now we are only making ideal pt2
    pt2["is_ideal"] = ak.ones_like(pt2.pls_pt, dtype=bool)

    return pt2
