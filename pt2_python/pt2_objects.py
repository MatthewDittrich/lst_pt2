import awkward as ak


def make_ls(events):
    ls = ak.zip ({
        "pt": events.ls_pt,
        "phi": events.ls_phi,
        "eta": events.ls_eta,
        "isFake": events.ls_isFake,
        "simIdx": events.ls_simIdx,
    })
    return ls

def make_pls(events):
    pls = ak.zip ({
        "pt": events.pls_pt,
        "phi": events.pls_phi,
        "eta", events.pls_eta,
        "nhit", events.pls_nhit,
        "isFake", events.pls_isFake,
        "simIdx", events.pls_simIdx,
    })
    return pls
