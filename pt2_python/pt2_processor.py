import awkward as ak

from tools import pt2_objects
from tools import pt2_methods

def pt2_processor(events):

    # Define the Sim Length
    #sim_length = ak.singletons(ak.num(events.sim_pt))
    sim_length = ak.num(events.sim_pt)

    # Initialize the LS and pLS objects
    ls = pt2_objects.get_ls(events)
    pls = pt2_objects.get_pls(events)

    # Attach the simIdx matching masks
    ls["matches_used_pv_pls"] = pt2_methods.ls_simIdx_matching(ls, pls, "used","pv", sim_length) 
    ls["matches_used_pu_pls"] = pt2_methods.ls_simIdx_matching(ls, pls, "used","pu", sim_length)
    ls["matches_unused_pv_pls"] = pt2_methods.ls_simIdx_matching(ls, pls, "unused","pv", sim_length)
    ls["matches_unused_pu_pls"] = pt2_methods.ls_simIdx_matching(ls, pls, "unused","pu", sim_length)
    pls["matches_used_pv_ls"] = pt2_methods.pls_simIdx_matching(pls, ls, "used","pv", sim_length) 
    pls["matches_used_pu_ls"] = pt2_methods.pls_simIdx_matching(pls, ls, "used","pu", sim_length)
    pls["matches_unused_pv_ls"] = pt2_methods.pls_simIdx_matching(pls, ls, "unused","pv", sim_length)
    pls["matches_unused_pu_ls"] = pt2_methods.pls_simIdx_matching(pls, ls, "unused","pu", sim_length)

    # simIdx mask combinations
    ls_matches_both_pv_pls = (ls.matches_used_pv_pls & ls.matches_unused_pv_pls)
    ls_matches_both_pu_pls = (ls.matches_used_pu_pls & ls.matches_unused_pu_pls)
    ls_matches_used_pv_pls = (ls.matches_used_pv_pls & (~ls_matches_both_pv_pls))
    ls_matches_unused_pv_pls = (ls.matches_unused_pv_pls & (~ls_matches_both_pv_pls))
    ls_matches_used_pu_pls = (ls.matches_used_pu_pls & (~ls_matches_both_pu_pls))
    ls_matches_unused_pu_pls = (ls.matches_unused_pu_pls & (~ls_matches_both_pu_pls))
    ls_matches_unused_pls = (ls.matches_unused_pv_pls | ls.matches_unused_pu_pls)
    ls_matches_used_pls = (ls.matches_used_pv_pls | ls.matches_used_pu_pls)
    ls_matches_none = (~(ls_matches_unused_pls | ls_matches_used_pls))
    pls_matches_both_pv_ls = (pls.matches_used_pv_ls & pls.matches_unused_pv_ls)
    pls_matches_both_pu_ls = (pls.matches_used_pu_ls & pls.matches_unused_pu_ls)
    pls_matches_used_pv_ls = (pls.matches_used_pv_ls & (~pls_matches_both_pv_ls))
    pls_matches_unused_pv_ls = (pls.matches_unused_pv_ls & (~pls_matches_both_pv_ls))
    pls_matches_used_pu_ls = (pls.matches_used_pu_ls & (~pls_matches_both_pu_ls))
    pls_matches_unused_pu_ls = (pls.matches_unused_pu_ls & (~pls_matches_both_pu_ls))
    pls_matches_unused_ls = (pls.matches_unused_pv_ls | pls.matches_unused_pu_ls)
    pls_matches_used_ls = (pls.matches_used_pv_ls | pls.matches_used_pu_ls)
    pls_matches_none = (~(pls_matches_unused_ls | pls_matches_used_ls))

    all_masks = {
        "ls": {
            "all_fake"                    : ls.isfake,
            "all_real"                    : ~ls.isfake,
            "used_fake"                   : (ls.isfake & ls.is_used),
            "used_real"                   : (~ls.isfake & ls.is_used),
            "unused_fake"                 : (ls.isfake & ~ls.is_used),
            "unused_real"                 : (~ls.isfake & ~ls.is_used),
            "simIdxmatches_none_pls"      : (ls_matches_none),
            "simIdxmatches_pv_used_pls"   : (ls_matches_used_pv_pls),
            "simIdxmatches_pv_unused_pls" : (ls_matches_unused_pv_pls),
            "simIdxmatches_pu_used_pls"   : (ls_matches_used_pu_pls),
            "simIdxmatches_pu_unused_pls" : (ls_matches_unused_pu_pls),
            "simIdxmatches_unused_pls"    : (ls_matches_unused_pls),
            "simIdxmatches_used_pls"      : (ls_matches_used_pls),
            "simIdxmatches_both_pv_pls"   : (ls_matches_both_pv_pls),
            "simIdxmatches_both_pu_pls"   : (ls_matches_both_pu_pls),
        },
        "pls": {
            "all_fake"                    : pls.isfake,
            "all_real"                    : ~pls.isfake,
            "used_fake"                   : (pls.isfake & pls.is_used),
            "used_real"                   : (~pls.isfake & pls.is_used),
            "unused_fake"                 : (pls.isfake & ~pls.is_used),
            "unused_real"                 : (~pls.isfake & ~pls.is_used),
            "simIdxmatches_none_ls"       : (pls_matches_none),
            "simIdxmatches_pv_used_ls"    : (pls_matches_used_pv_ls),
            "simIdxmatches_pv_unused_ls"  : (pls_matches_unused_pv_ls),
            "simIdxmatches_pu_used_ls"    : (pls_matches_used_pu_ls),
            "simIdxmatches_pu_unused_ls"  : (pls_matches_unused_pu_ls),
            "simIdxmatches_unused_ls"     : (pls_matches_unused_ls),
            "simIdxmatches_used_ls"       : (pls_matches_used_ls),
            "simIdxmatches_both_pv_ls"    : (pls_matches_both_pv_ls),
            "simIdxmatches_both_pu_ls"    : (pls_matches_both_pu_ls),
        }
    }

    ###############################
    ### Making the Plot Objects ###
    ###############################

    # Define the empty plot objects
    plot_objects = {}

    # Map to loop through objects
    objects = {
        "ls"  : ls,
        "pls" : pls,
    }

    # Commonly plotted values
    common_values = ["pt","eta","phi"]    

    for obj_key, obj in objects.items():
        for cat, mask in all_masks[obj_key].items():
            if cat.startswith("simIdxmatches"):
                # simIdx matching is done for used/real and unused/real
                for val in common_values:
                    used_real_key = f"{obj_key}_used_real_{cat}_{val}"
                    used_real_obj = getattr(obj,val)[mask & all_masks[obj_key]["used_real"]]
                    unused_real_key = f"{obj_key}_unused_real_{cat}_{val}"
                    unused_real_obj = getattr(obj,val)[mask & all_masks[obj_key]["unused_real"]]

                    plot_objects[used_real_key] = used_real_obj
                    plot_objects[unused_real_key] = unused_real_obj
            else:
                # Right now we are just doing common values
                for val in common_values:
                    plot_objects_key = f"{obj_key}_{cat}_{val}"
                    plot_objects_array = getattr(obj,val)[mask]                                    

                    plot_objects[plot_objects_key] = plot_objects_array

    return plot_objects
