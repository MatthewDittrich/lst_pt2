import uproot
import awkward as ak
import argparse

import pt2_objects

def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file_path", help = "The path to the root file")
    parser.add_argument("-o", "--ouput-path", default = "pt2_output", help = "The path the output will be saved to")
    args = parser.args()

    #Load the root file and grab the events
    rfile = uproot.open(args.root_file_path)
    if 'tree' not in rfile:
        raise Exception("Tree 'tree' was not found in the root file! Please ensure you are using an output from LST OD.")
    else:
        tree = rfile["tree"]

    # Grab the events
    events = tree.arrays(entry_start=0, entry_stop=(tree.num_entries), library="ak")

    # Make some wanted objects
    ls = pt2_objects.make_ls(events)
    pls = pt2_objects.make_pls(events)

    # Make the all important used/unused masks
    used_ls_mask, used_pls_mask = pt2_methods.ls_pls_used_mask() 
