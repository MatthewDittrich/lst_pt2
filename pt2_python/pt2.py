import uproot
import awkward as ak
import argparse
import multiprocessing
import os
import concurrent.futures  # Updated to use concurrent.futures

from tools import pt2_methods
from tools import pt2_plotting
from tools.GatorMaker import print_gator


branches_to_load = [
    "sim_pt",
    "ls_pt", "ls_eta", "ls_phi", "ls_isFake", "ls_simIdx",
    "pls_pt", "pls_eta", "pls_phi", "pls_isFake", "pls_simIdx",
    "tc_pt5Idx", "tc_t5Idx", "tc_pt3Idx","tc_plsIdx",
    "t3_lsIdx0","t3_lsIdx1",
    "t5_t3Idx0","t5_t3Idx1",
    "pt3_plsIdx","pt3_t3Idx",
    "pt5_plsIdx","pt5_t5Idx",
]



def pt2_processor(events):

    # Make the used masks
    used_ls_mask, used_pls_mask = pt2_methods.ls_pls_used_masks(events)

    # Make the fake masks
    fake_ls_mask = (events.ls_isFake == 1)
    fake_pls_mask = (events.pls_isFake == 1)

    # Define some of the final masks
    # LS
    fake_used_ls_mask = (fake_ls_mask & used_ls_mask)
    real_used_ls_mask = ((~fake_ls_mask) & used_ls_mask)
    fake_unused_ls_mask = (fake_ls_mask & (~used_ls_mask))
    real_unused_ls_mask = ((~fake_ls_mask) & (~used_ls_mask))
    # pLS
    fake_used_pls_mask = (fake_pls_mask & used_pls_mask)
    real_used_pls_mask = ((~fake_pls_mask) & used_pls_mask)
    fake_unused_pls_mask = (fake_pls_mask & (~used_pls_mask))
    real_unused_pls_mask = ((~fake_pls_mask) & (~used_pls_mask))

    # Make masks for simIdx matching
    ls_matches_usedrealpls_mask, \
    ls_matches_unusedrealpls_mask = pt2_methods.ls_pls_matching_simIdx(
        events, 
        used_ls_mask, 
        fake_ls_mask, 
        used_pls_mask, 
        fake_pls_mask, 
        "ls"
    ) 

    pls_matches_usedrealls_mask, \
    pls_matches_unusedrealls_mask = pt2_methods.ls_pls_matching_simIdx(
        events, 
        used_ls_mask, 
        fake_ls_mask, 
        used_pls_mask, 
        fake_pls_mask, 
        "pls"
    ) 


    # Final objects dictionary
    plot_objects = {
        "ls_all_isfake_pt": events.ls_pt[fake_ls_mask],
        "ls_all_isfake_phi": events.ls_phi[fake_ls_mask],
        "ls_all_isfake_eta": events.ls_eta[fake_ls_mask],
        "ls_all_real_pt": events.ls_pt[~fake_ls_mask],
        "ls_all_real_phi": events.ls_phi[~fake_ls_mask],
        "ls_all_real_eta": events.ls_eta[~fake_ls_mask],
        "ls_used_isfake_pt": events.ls_pt[fake_used_ls_mask],
        "ls_used_isfake_eta": events.ls_eta[fake_used_ls_mask],
        "ls_used_isfake_phi": events.ls_phi[fake_used_ls_mask],
        "ls_used_real_pt": events.ls_pt[real_used_ls_mask],
        "ls_used_real_eta": events.ls_eta[real_used_ls_mask],
        "ls_used_real_phi": events.ls_phi[real_used_ls_mask],
        "ls_unused_isfake_pt": events.ls_pt[fake_unused_ls_mask],
        "ls_unused_isfake_eta": events.ls_eta[fake_unused_ls_mask],
        "ls_unused_isfake_phi": events.ls_phi[fake_unused_ls_mask],
        "ls_unused_real_pt": events.ls_pt[real_unused_ls_mask],
        "ls_unused_real_eta": events.ls_eta[real_unused_ls_mask],
        "ls_unused_real_phi": events.ls_phi[real_unused_ls_mask],
        "ls_used_real_simIdxmacthes_pls_used_real_pt": events.ls_pt[real_used_ls_mask & ls_matches_usedrealpls_mask],
        "ls_used_real_simIdxmacthes_pls_used_real_eta": events.ls_eta[real_used_ls_mask & ls_matches_usedrealpls_mask],
        "ls_used_real_simIdxmacthes_pls_used_real_phi": events.ls_phi[real_used_ls_mask & ls_matches_usedrealpls_mask],
        "ls_used_real_simIdxmacthes_pls_unused_real_pt": events.ls_pt[real_used_ls_mask & ls_matches_unusedrealpls_mask],
        "ls_used_real_simIdxmacthes_pls_unused_real_eta": events.ls_eta[real_used_ls_mask & ls_matches_unusedrealpls_mask],
        "ls_used_real_simIdxmacthes_pls_unused_real_phi": events.ls_phi[real_used_ls_mask & ls_matches_unusedrealpls_mask],

        "ls_used_real_simIdxmacthes_none_pt": events.ls_pt[real_used_ls_mask & ~ls_matches_unusedrealpls_mask & ~ls_matches_usedrealpls_mask],
        "ls_used_real_simIdxmacthes_none_eta": events.ls_eta[real_used_ls_mask & ~ls_matches_unusedrealpls_mask & ~ls_matches_usedrealpls_mask],
        "ls_used_real_simIdxmacthes_none_phi": events.ls_phi[real_used_ls_mask & ~ls_matches_unusedrealpls_mask & ~ls_matches_usedrealpls_mask],

        "ls_unused_real_simIdxmacthes_pls_used_real_pt": events.ls_pt[real_unused_ls_mask & ls_matches_usedrealpls_mask],
        "ls_unused_real_simIdxmacthes_pls_used_real_eta": events.ls_eta[real_unused_ls_mask & ls_matches_usedrealpls_mask],
        "ls_unused_real_simIdxmacthes_pls_used_real_phi": events.ls_phi[real_unused_ls_mask & ls_matches_usedrealpls_mask],
        "ls_unused_real_simIdxmacthes_pls_unused_real_pt": events.ls_pt[real_unused_ls_mask & ls_matches_unusedrealpls_mask],
        "ls_unused_real_simIdxmacthes_pls_unused_real_eta": events.ls_eta[real_unused_ls_mask & ls_matches_unusedrealpls_mask],
        "ls_unused_real_simIdxmacthes_pls_unused_real_phi": events.ls_phi[real_unused_ls_mask & ls_matches_unusedrealpls_mask],

        "ls_unused_real_simIdxmacthes_none_pt": events.ls_pt[real_unused_ls_mask & (~ls_matches_unusedrealpls_mask) & (~ls_matches_usedrealpls_mask)],
        "ls_unused_real_simIdxmacthes_none_eta": events.ls_eta[real_unused_ls_mask & (~ls_matches_unusedrealpls_mask) & (~ls_matches_usedrealpls_mask)],
        "ls_unused_real_simIdxmacthes_none_phi": events.ls_phi[real_unused_ls_mask & (~ls_matches_unusedrealpls_mask) & (~ls_matches_usedrealpls_mask)],

        "pls_all_isfake_pt": events.pls_pt[fake_pls_mask],
        "pls_all_isfake_phi": events.pls_phi[fake_pls_mask],
        "pls_all_isfake_eta": events.pls_eta[fake_pls_mask],
        "pls_all_real_pt": events.pls_pt[~fake_pls_mask],
        "pls_all_real_phi": events.pls_phi[~fake_pls_mask],
        "pls_all_real_eta": events.pls_eta[~fake_pls_mask],
        "pls_used_isfake_pt": events.pls_pt[fake_used_pls_mask],
        "pls_used_isfake_eta": events.pls_eta[fake_used_pls_mask],
        "pls_used_isfake_phi": events.pls_phi[fake_used_pls_mask],
        "pls_used_real_pt": events.pls_pt[real_used_pls_mask],
        "pls_used_real_eta": events.pls_eta[real_used_pls_mask],
        "pls_used_real_phi": events.pls_phi[real_used_pls_mask],
        "pls_unused_isfake_pt": events.pls_pt[fake_unused_pls_mask],
        "pls_unused_isfake_eta": events.pls_eta[fake_unused_pls_mask],
        "pls_unused_isfake_phi": events.pls_phi[fake_unused_pls_mask],
        "pls_unused_real_pt": events.pls_pt[real_unused_pls_mask],
        "pls_unused_real_eta": events.pls_eta[real_unused_pls_mask],
        "pls_unused_real_phi": events.pls_phi[real_unused_pls_mask],
        "pls_used_real_simIdxmacthes_ls_used_real_pt": events.pls_pt[real_used_pls_mask & pls_matches_usedrealls_mask],
        "pls_used_real_simIdxmacthes_ls_used_real_eta": events.pls_eta[real_used_pls_mask & pls_matches_usedrealls_mask],
        "pls_used_real_simIdxmacthes_ls_used_real_phi": events.pls_phi[real_used_pls_mask & pls_matches_usedrealls_mask],
        "pls_used_real_simIdxmacthes_ls_unused_real_pt": events.pls_pt[real_used_pls_mask & pls_matches_unusedrealls_mask],
        "pls_used_real_simIdxmacthes_ls_unused_real_eta": events.pls_eta[real_used_pls_mask & pls_matches_unusedrealls_mask],
        "pls_used_real_simIdxmacthes_ls_unused_real_phi": events.pls_phi[real_used_pls_mask & pls_matches_unusedrealls_mask],
        "pls_used_real_simIdxmacthes_none_pt": events.pls_pt[real_used_pls_mask & (~pls_matches_unusedrealls_mask) & (~pls_matches_usedrealls_mask)],
        "pls_used_real_simIdxmacthes_none_eta": events.pls_eta[real_used_pls_mask & (~pls_matches_unusedrealls_mask) & (~pls_matches_usedrealls_mask)],
        "pls_used_real_simIdxmacthes_none_phi": events.pls_phi[real_used_pls_mask & (~pls_matches_unusedrealls_mask) & (~pls_matches_usedrealls_mask)],
        "pls_unused_real_simIdxmacthes_ls_used_real_pt": events.pls_pt[real_unused_pls_mask & pls_matches_usedrealls_mask],
        "pls_unused_real_simIdxmacthes_ls_used_real_eta": events.pls_eta[real_unused_pls_mask & pls_matches_usedrealls_mask],
        "pls_unused_real_simIdxmacthes_ls_used_real_phi": events.pls_phi[real_unused_pls_mask & pls_matches_usedrealls_mask],
        "pls_unused_real_simIdxmacthes_ls_unused_real_pt": events.pls_pt[real_unused_pls_mask & pls_matches_unusedrealls_mask],
        "pls_unused_real_simIdxmacthes_ls_unused_real_eta": events.pls_eta[real_unused_pls_mask & pls_matches_unusedrealls_mask],
        "pls_unused_real_simIdxmacthes_ls_unused_real_phi": events.pls_phi[real_unused_pls_mask & pls_matches_unusedrealls_mask],
        "pls_unused_real_simIdxmacthes_none_pt": events.pls_pt[real_unused_pls_mask & (~pls_matches_unusedrealls_mask) & (~pls_matches_usedrealls_mask)],
        "pls_unused_real_simIdxmacthes_none_eta": events.pls_eta[real_unused_pls_mask & (~pls_matches_unusedrealls_mask) & (~pls_matches_usedrealls_mask)],
        "pls_unused_real_simIdxmacthes_none_phi": events.pls_phi[real_unused_pls_mask & (~pls_matches_unusedrealls_mask) & (~pls_matches_usedrealls_mask)],
    }

    return plot_objects


def merge_results(results_list):
    merged = {}
    keys = results_list[0].keys()
    for key in keys:
        arrays_to_concat = [r[key] for r in results_list]
        merged[key] = ak.concatenate(arrays_to_concat)
    return merged


def main():

    print_gator()

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file_path", help="The path to the root file")
    parser.add_argument("-o", "--output-path", default="pt2_output", help="The path the output will be saved to")
    parser.add_argument("-n", "--number-workers", type=int, default=max(multiprocessing.cpu_count() - 2, 1),help="The number of workers")
    parser.add_argument("-c", "--chunk-size", type=int, default=5, help="Chunk Size")
    args = parser.parse_args()

    print("Attempting to get events...")

    # Set up ROOT chunk iterator
    events_iter = uproot.iterate(f"{args.root_file_path}:tree", step_size=args.chunk_size, filter_name=branches_to_load, library="ak")
    print("Events will be chunked into groups of", args.chunk_size)
    chunk_list = list(events_iter)
    print(f"Total chunks: {len(chunk_list)}")

    # Use concurrent.futures for parallel processing
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.number_workers) as executor:
        results_list = list(executor.map(pt2_processor, chunk_list))

    # Merge the results
    merged_results = merge_results(results_list)

    # Now pass merged_results to plotting function
    pt2_plotting.make_plots(merged_results, args.output_path)
    print("Plots Finished!")


if __name__ == "__main__":
    main()
