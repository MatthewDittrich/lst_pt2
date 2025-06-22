import awkward as ak

from tools import pt2_methods


class LST_OBJECT:
    
    def __init__(self, data: dict):
        self._data = data
    
    def __getattr__(self, key):
        try:
            return self._data[key]
        except KeyError:
            raise AttributeError(f"'LST_OBJECT' has no attribute '{key}'")

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        if not self._data:
            raise ValueError("Cannot add a field to an empty LST_OBJECT.")
        # Reference key to compare shapes against
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
