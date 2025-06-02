import awkward as ak

def unique_per_event(array_of_arrays):
    return ak.Array([list(set(x)) for x in array_of_arrays])


def ls_pls_matching_simIdx(events, ls_isused_mask, ls_isfake_mask, pls_isused_mask, pls_isfake_mask, key):
    # Goal of this function is to create masks
    # First we will create ls masks for having a sim matched fake/real and used/unused pls
    # Then we will create a pls mask for having a sim matched fake/real and used/unused ls

    if key not in ["ls","pls"]:
        raise Exception(f"Key {key} is not recognized!")
    elif (key == "ls"):
        # Define the LS simIdx
        ls_simIdx = events.ls_simIdx
        # Get the filtered pLS simIdx
        used_real_pls = events.pls_simIdx[pls_isused_mask & (~pls_isfake_mask) & (events.pls_simIdx != -999) & (events.pls_simIdx < ak.num(events.sim_pt))]
        unused_real_pls = events.pls_simIdx[(~pls_isused_mask) & (~pls_isfake_mask) & (events.pls_simIdx != -999) & (events.pls_simIdx < ak.num(events.sim_pt))]
        # Construct the masks
        ls_matches_usedrealpls_mask = []
        ls_matches_unusedrealpls_mask = []
        # Loop over the lists
        # Not broadcasting due to memory issues
        for ls_list, pls_usedreal_list, pls_unusedreal_list in zip(ls_simIdx, used_real_pls, unused_real_pls):
            pls_usedreal_set = set(pls_usedreal_list)
            pls_unusedreal_set = set(pls_unusedreal_list) 
            ls_matches_usedrealpls_mask.append([a in pls_usedreal_set for a in ls_list])
            ls_matches_unusedrealpls_mask.append([a in pls_unusedreal_set for a in ls_list])

        # Return the masks
        return ak.Array(ls_matches_usedrealpls_mask), ak.Array(ls_matches_unusedrealpls_mask)
 
    elif (key == "pls"):
        # Define the pLS simIdx
        pls_simIdx = events.pls_simIdx
        # Get the filtered LS simIdx
        used_real_ls = events.ls_simIdx[ls_isused_mask & (~ls_isfake_mask) & (events.ls_simIdx != -999) & (events.ls_simIdx < ak.num(events.sim_pt))]
        unused_real_ls = events.ls_simIdx[(~ls_isused_mask) & (~ls_isfake_mask) & (events.ls_simIdx != -999) & (events.ls_simIdx < ak.num(events.sim_pt))]
        # Construct the masks
        pls_matches_usedrealls_mask = []
        pls_matches_unusedrealls_mask = []
        # Loop over the lists
        # Not broadcasting due to memory issues
        for pls_list, ls_usedreal_list, ls_unusedreal_list in zip(pls_simIdx, used_real_ls, unused_real_ls):
            ls_usedreal_set = set(ls_usedreal_list)
            ls_unusedreal_set = set(ls_unusedreal_list) 
            pls_matches_usedrealls_mask.append([a in ls_usedreal_set for a in pls_list])
            pls_matches_unusedrealls_mask.append([a in ls_unusedreal_set for a in pls_list])

        # Return the masks
        return ak.Array(pls_matches_usedrealls_mask), ak.Array(pls_matches_unusedrealls_mask)


def ls_pls_used_masks(events):

    # Grab Idx for each TC type
    tc_pls_valid = (events.tc_plsIdx != -999)
    tc_pt3_valid = (events.tc_pt3Idx != -999)
    tc_pt5_valid = (events.tc_pt5Idx != -999)
    tc_t5_valid  = (events.tc_t5Idx != -999)

    # Gather pLS indices used directly in TC
    pls_direct = events.tc_plsIdx[tc_pls_valid]

    # pT3: pLS and LS
    pt3_idx = events.tc_pt3Idx[tc_pt3_valid]
    pls_from_pt3 = events.pt3_plsIdx[pt3_idx]
    t3_from_pt3 = events.pt3_t3Idx[pt3_idx]
    ls0_from_pt3 = events.t3_lsIdx0[t3_from_pt3]
    ls1_from_pt3 = events.t3_lsIdx1[t3_from_pt3]

    # pT5: pLS and LS
    pt5_idx = events.tc_pt5Idx[tc_pt5_valid]
    pls_from_pt5 = events.pt5_plsIdx[pt5_idx]
    t5_from_pt5 = events.pt5_t5Idx[pt5_idx]
    t3_0_from_pt5 = events.t5_t3Idx0[t5_from_pt5]
    t3_1_from_pt5 = events.t5_t3Idx1[t5_from_pt5]
    ls0_from_pt5 = events.t3_lsIdx0[t3_0_from_pt5]
    ls1_from_pt5 = events.t3_lsIdx1[t3_0_from_pt5]
    ls2_from_pt5 = events.t3_lsIdx0[t3_1_from_pt5]
    ls3_from_pt5 = events.t3_lsIdx1[t3_1_from_pt5]

    # T5: LS
    t5_idx = events.tc_t5Idx[tc_t5_valid]
    t3_0_from_t5 = events.t5_t3Idx0[t5_idx]
    t3_1_from_t5 = events.t5_t3Idx1[t5_idx]
    ls0_from_t5 = events.t3_lsIdx0[t3_0_from_t5]
    ls1_from_t5 = events.t3_lsIdx1[t3_0_from_t5]
    ls2_from_t5 = events.t3_lsIdx0[t3_1_from_t5]
    ls3_from_t5 = events.t3_lsIdx1[t3_1_from_t5]

    # Concatenate all the used pLS and LS
    all_used_pls = ak.concatenate([pls_direct, pls_from_pt3, pls_from_pt5], axis=1)
    all_used_ls = ak.concatenate([
        ls0_from_pt3, ls1_from_pt3,
        ls0_from_pt5, ls1_from_pt5, ls2_from_pt5, ls3_from_pt5,
        ls0_from_t5, ls1_from_t5, ls2_from_t5, ls3_from_t5,
    ], axis=1)

    # Remove duplicates
    # TODO: Is this needed?
    unique_used_pls = unique_per_event(all_used_pls)
    unique_used_ls  = unique_per_event(all_used_ls)

    # Full index arrays of LS and pLS
    full_pls_idx = ak.local_index(events.pls_pt)
    full_ls_idx  = ak.local_index(events.ls_pt)

    # Create masks for used LS and pLS
    pls_mask_list = []
    for full, used in zip(full_pls_idx, unique_used_pls):
        used_set = set(used)
        pls_mask_list.append([idx in used_set for idx in full])
    pls_mask = ak.Array(pls_mask_list)

    ls_mask_list = []
    for full, used in zip(full_ls_idx, unique_used_ls):
        used_set = set(used)
        ls_mask_list.append([idx in used_set for idx in full])
    ls_mask = ak.Array(ls_mask_list)

    return ls_mask, pls_mask
