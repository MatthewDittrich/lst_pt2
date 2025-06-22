import awkward as ak



def unique_per_event(array_of_arrays):
    return ak.Array([list(set(x)) for x in array_of_arrays])



def ls_simIdx_matching(ls, pls, used_key, hard_key, sim_length):
    if hard_key not in ["pv", "pu"]:
        raise Exception(f"Key {hard_key} is not recognized!")
    if used_key not in ["used","unused"]:
        raise Exception(f"Key {used_key} is not recognized!")
    # Define the simIdx and needed masks
    ls_simIdx = ls["simIdx"]
    pls_simIdx = pls["simIdx"]
    real_pls_mask = (~pls["isfake"])
    pls_used_mask = pls["is_used"]
    # Base mask for filtering pls_simidx
    base_mask = (real_pls_mask & (pls_simIdx != -999))
    # Filters are dependent on keys
    if (hard_key == "pv") and (used_key == "used"):
        base_mask = base_mask & pls_used_mask & (pls_simIdx < sim_length)
    elif (hard_key == "pv") and (used_key == "unused"):
        base_mask = base_mask & (~pls_used_mask) & (pls_simIdx < sim_length)
    elif (hard_key == "pu") and (used_key == "used"):
        base_mask = base_mask & pls_used_mask & (pls_simIdx >= sim_length)
    elif (hard_key == "pu") and (used_key == "unused"):
        base_mask = base_mask & (~pls_used_mask) & (pls_simIdx >= sim_length)
    pls_simIdx = pls_simIdx[base_mask]
    # Construct the mask
    ls_mask = []
    for ls_list, pls_list in zip(ls_simIdx, pls_simIdx):
        pls_set = set(pls_list) 
        ls_mask.append([a in pls_set for a in ls_list])
    # Return the mask
    return ak.Array(ls_mask)



def pls_simIdx_matching(pls, ls, used_key, hard_key, sim_length):
    if hard_key not in ["pv", "pu"]:
        raise Exception(f"Key {hard_key} is not recognized!")
    if used_key not in ["used","unused"]:
        raise Exception(f"Key {used_key} is not recognized!")
    # Define the simIdx and needed masks
    ls_simIdx = ls["simIdx"]
    pls_simIdx = pls["simIdx"]
    real_ls_mask = (~ls["isfake"])
    ls_used_mask = ls["is_used"]
    # Base mask for filtering pls_simidx
    base_mask = (real_ls_mask & (ls_simIdx != -999))
    # Filters are dependent on keys
    if (hard_key == "pv") and (used_key == "used"):
        base_mask = base_mask & ls_used_mask & (ls_simIdx < sim_length)
    elif (hard_key == "pv") and (used_key == "unused"):
        base_mask = base_mask & (~ls_used_mask) & (ls_simIdx < sim_length)
    elif (hard_key == "pu") and (used_key == "used"):
        base_mask = base_mask & ls_used_mask & (ls_simIdx >= sim_length)
    elif (hard_key == "pu") and (used_key == "unused"):
        base_mask = base_mask & (~ls_used_mask) & (ls_simIdx >= sim_length)
    ls_simIdx = ls_simIdx[base_mask]
    # Construct the mask
    pls_mask = []
    for pls_list, ls_list in zip(pls_simIdx, ls_simIdx):
        ls_set = set(ls_list) 
        pls_mask.append([a in ls_set for a in pls_list])
    # Return the mask
    return ak.Array(pls_mask)



def ls_masks(events, key):

    # Full index arrays for LS
    full_ls_idx  = ak.local_index(events.ls_pt)

    if key not in ["pt5", "t5", "pt3"]:
        raise Exception(f"LS Key {key} is not a recognized option!")

    elif (key == "pt5"):
        # Grab the pt5 Idxs used in TCs
        pt5_idx = events.tc_pt5Idx[events.tc_pt5Idx != -999]
        # breakdown pt5 into 4 LS
        t5_from_pt5 = events.pt5_t5Idx[pt5_idx]
        t3_0_from_pt5 = events.t5_t3Idx0[t5_from_pt5]
        t3_1_from_pt5 = events.t5_t3Idx1[t5_from_pt5]
        ls0_from_pt5 = events.t3_lsIdx0[t3_0_from_pt5]
        ls1_from_pt5 = events.t3_lsIdx1[t3_0_from_pt5]
        ls2_from_pt5 = events.t3_lsIdx0[t3_1_from_pt5]
        ls3_from_pt5 = events.t3_lsIdx1[t3_1_from_pt5]
        # Concatenate 
        all_pt5_ls = ak.concatenate([
            ls0_from_pt5, ls1_from_pt5, ls2_from_pt5, ls3_from_pt5,
        ], axis=1)
        # Remove duplicates
        # TODO: Is this needed?
        unique_pt5_ls = unique_per_event(all_pt5_ls)
        # Create the Mask
        ls_from_pt5_mask_list = []
        for full, pt5 in zip(full_ls_idx, unique_pt5_ls):
            pt5_set = set(pt5)
            ls_from_pt5_mask_list.append([idx in pt5_set for idx in full])
        ls_from_pt5_mask = ak.Array(ls_from_pt5_mask_list)
        
        return ls_from_pt5_mask

    elif (key == "t5"):
        # Grab the t5 Idxs used in TCs
        t5_idx = events.tc_t5Idx[events.tc_t5Idx != -999]
        # breakdown t5 into 4 LS
        t3_0_from_t5 = events.t5_t3Idx0[t5_idx]
        t3_1_from_t5 = events.t5_t3Idx1[t5_idx]
        ls0_from_t5 = events.t3_lsIdx0[t3_0_from_t5]
        ls1_from_t5 = events.t3_lsIdx1[t3_0_from_t5]
        ls2_from_t5 = events.t3_lsIdx0[t3_1_from_t5]
        ls3_from_t5 = events.t3_lsIdx1[t3_1_from_t5]
        # Concatenate 
        all_t5_ls = ak.concatenate([
            ls0_from_t5, ls1_from_t5, ls2_from_t5, ls3_from_t5,
        ], axis=1)
        # Remove duplicates
        # TODO: Is this needed?
        unique_t5_ls = unique_per_event(all_t5_ls)
        # Create the Mask
        ls_from_t5_mask_list = []
        for full, t5 in zip(full_ls_idx, unique_t5_ls):
            t5_set = set(t5)
            ls_from_t5_mask_list.append([idx in t5_set for idx in full])
        ls_from_t5_mask = ak.Array(ls_from_t5_mask_list)

        return ls_from_t5_mask

    elif (key == "pt3"):
        # Grab the pt3 Idxs used in TCs
        pt3_idx = events.tc_pt3Idx[events.tc_pt3Idx != -999]
        # breakdown pt3 into 2 LS
        t3_from_pt3 = events.pt3_t3Idx[pt3_idx]
        ls0_from_pt3 = events.t3_lsIdx0[t3_from_pt3]
        ls1_from_pt3 = events.t3_lsIdx1[t3_from_pt3]
        # Concatenate 
        all_pt3_ls = ak.concatenate([
            ls0_from_pt3, ls1_from_pt3,
        ], axis=1)
        # Remove duplicates
        # TODO: Is this needed?
        unique_pt3_ls = unique_per_event(all_pt3_ls)
        # Create the Mask
        ls_from_pt3_mask_list = []
        for full, pt3 in zip(full_ls_idx, unique_pt3_ls):
            pt3_set = set(pt3)
            ls_from_pt3_mask_list.append([idx in pt3_set for idx in full])
        ls_from_pt3_mask = ak.Array(ls_from_pt3_mask_list)

        return ls_from_pt3_mask



def pls_masks(events, key):

    # Full index arrays for pLS
    full_pls_idx  = ak.local_index(events.pls_pt)

    if key not in ["pt5", "pLS", "pt3"]:
        raise Exception(f"Key {key} is not a recognized option!")

    elif (key == "pt5"):
        # Grab the pt5 Idxs used in TCs
        pt5_idx = events.tc_pt5Idx[events.tc_pt5Idx != -999]
        # Grab the pLS from pt5
        pls_from_pt5 = events.pt5_plsIdx[pt5_idx]
        # Remove duplicates
        # TODO: Is this needed?
        unique_pt5_pls = unique_per_event(pls_from_pt5)
        # Create the Mask
        pls_from_pt5_mask_list = []
        for full, pt5 in zip(full_pls_idx, unique_pt5_pls):
            pt5_set = set(pt5)
            pls_from_pt5_mask_list.append([idx in pt5_set for idx in full])
        pls_from_pt5_mask = ak.Array(pls_from_pt5_mask_list)
        
        return pls_from_pt5_mask

    elif (key == "pt3"):
        # Grab the pt3 Idxs used in TCs
        pt3_idx = events.tc_pt3Idx[events.tc_pt3Idx != -999]
        # Grab the pLS from pt3
        pls_from_pt3 = events.pt3_plsIdx[pt3_idx]
        # Remove duplicates
        # TODO: Is this needed?
        unique_pt3_pls = unique_per_event(pls_from_pt3)
        # Create the Mask
        pls_from_pt3_mask_list = []
        for full, pt3 in zip(full_pls_idx, unique_pt3_pls):
            pt3_set = set(pt3)
            pls_from_pt3_mask_list.append([idx in pt3_set for idx in full])
        pls_from_pt3_mask = ak.Array(pls_from_pt3_mask_list)
        
        return pls_from_pt3_mask

    elif (key == "pLS"):
        # Grab the pLS Idxs used in TCs
        pls_idx = events.tc_plsIdx[events.tc_plsIdx != -999]
        # Remove duplicates
        # TODO: Is this needed?
        unique_tc_pls = unique_per_event(pls_idx)
        # Create the Mask
        pls_from_tc_mask_list = []
        for full, pls in zip(full_pls_idx, unique_tc_pls):
            pls_set = set(pls)
            pls_from_tc_mask_list.append([idx in pls_set for idx in full])
        pls_from_tc_mask = ak.Array(pls_from_tc_mask_list)
        
        return pls_from_tc_mask
