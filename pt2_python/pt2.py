import uproot
import awkward as ak
import argparse
import multiprocessing
import concurrent.futures

from tools.branches import branches_to_load
from tools.plotting import pt2_plotting
from tools.GatorMaker import print_gator
from pt2_processor import pt2_processor

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
    parser.add_argument("-p", "--plot", action="store_true", help="Make plots")
    parser.add_argument("-y", "--yields", action="store_true", help="Make yield table")
    parser.add_argument("-r", "--root-file", action="store_true", help="Outputs root file with additional branches from the original")
    parser.add_argument("-n", "--number-workers", type=int, default=max(multiprocessing.cpu_count() - 2, 1), help="The number of workers")
    parser.add_argument("-c", "--chunk-size", type=int, default=5, help="Chunk Size")
    parser.add_argument("-d", "--debug", action="store_true", help="Run in debug mode (process only first event of first chunk)")
    args = parser.parse_args()

    # The code should be set up to output something
    if not (args.plot or args.root_file or args.yields):
        raise Exception("Code will not output anything! Please ensure some form of output.")

    print("Attempting to get events...")

    # Set up ROOT chunk iterator
    events_iter = uproot.iterate(f"{args.root_file_path}:tree", step_size=args.chunk_size, filter_name=branches_to_load, library="ak")
    print("Events will be chunked into groups of", args.chunk_size)
    chunk_list = list(events_iter)
    print(f"Total chunks: {len(chunk_list)}")

    if args.debug:
        print("DEBUG MODE ENABLED: Processing only 1 event from the first chunk.")
        first_chunk = chunk_list[0][:1]
        results_list = [pt2_processor(first_chunk)]
    else:
        # Use concurrent.futures for parallel processing
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.number_workers) as executor:
            results_list = list(executor.map(pt2_processor, chunk_list))

    # Merge the results
    merged_results = merge_results(results_list)

    if args.yields:
        # TODO Set this up
        raise Exception("Not set up yet!")

    if args.root_file:
        # TODO Set this up
        raise Exception("Not set up yet!")

    if args.plot:
        # Pass Merged Results to Plotting Function
        pt2_plotting.make_plots(merged_results, args.output_path)
        print("Plots Finished!")

if __name__ == "__main__":
    main()
