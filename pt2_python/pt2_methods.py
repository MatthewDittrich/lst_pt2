import awkward as ak

def ls_pls_used_mask (events):
# This method will return two masks based on whether a LS or pLS was used in a TC
# Each track candidate is broken down into LS and pLS

    # Define empty awkward arrays for used 
    used_ls = ak.Array([[]] * len(events))
    used_pls = ak.Array([[]] * len(events)) 

    # TC: pLS
    used_pls_idx = events.tc_plsIdx[events.tc_plsIdx != -999]
    used_pls = ak.concatenate([used_pls, used_pls_idx], axis = 1)

    # TC: pT3
    used_pt3_idx = events.tc_pt3Idx[events.tc_pt3Idx != -999]
    # First grab the pLS
    pls_from_pt3 = events.pt3_plsIdx[used_pt3_idx]
    used_pls = ak.concatenate([used_pls, pls_from_pt3], axis = 1)
    # Now grab the two LS
    t3_from_pt3 = events.pt3_t3Idx[used_pt3_idx]
    ls0_from_pt3 = events.t3_lsIdx0[t3_from_pt3]
    ls1_from_pt3 = events.t3_lsIdx1[t3_from_pt3]
    used_ls = ak.concatenate([used_ls, ls0_from_pt3, ls1_from_pt3], axis = 1)

    # TC: pT5
    used_pt5_idx = events.tc_pt5Idx[events.tc_pt5Idx != -999]
    # First grab the pLS
    pls_from_pt5 = events.pt5_plsIdx[used_pt5_idx]
    used_pls = ak.concatenate([used_pls, pls_from_pt5], axis = 1)
    # Grab the t5 from the pt5
    t5_from_pt5 = events.pt5_t5Idx[used_pt5_idx]
    # Now grab the two T3
    t30_from_pt5 = events.t5_t3Idx0[t5_from_pt5]
    t31_from_pt5 = events.t5_t3Idx1[t5_from_pt5]
    # Now grab the four LS
    ls0_from_pt5 = events.t3_lsIdx0[t30_from_pt5]
    ls1_from_pt5 = events.t3_lsIdx1[t30_from_pt5]
    ls2_from_pt5 = events.t3_lsIdx0[t31_from_pt5]
    ls3_from_pt5 = events.t3_lsIdx1[t31_from_pt5]
    used_ls = ak.concatenate([used_ls, ls0_from_pt5, ls1_from_pt5, ls2_from_pt5, ls3_from_pt5], axis = 1)

    # TC: T5
    used_t5_idx = events.tc_t5Idx[events.tc_t5Idx != -999]
    # Grab the two T3
    t30_from_t5 = events.t5_t3Idx0[used_t5_idx]
    t31_from_t5 = events.t5_t3Idx1[used_t5_idx]
    # Grab the four LS
    ls0_from_t5 = events.t3_lsIdx0[t30_from_t5]
    ls1_from_t5 = events.t3_lsIdx1[t30_from_t5]
    ls2_from_t5 = events.t3_lsIdx0[t31_from_t5]
    ls3_from_t5 = events.t3_lsIdx1[t31_from_t5]
    used_ls = ak.concatenate([used_ls, ls0_from_t5, ls1_from_t5, ls2_from_t5, ls3_from_t5], axis = 1)

    #

def unique_inner(ak_array):
    unique_lists = [np.unique(ak.to_numpy(sublist)) for sublist in ak_array]
    return ak.Array(unique_lists)
