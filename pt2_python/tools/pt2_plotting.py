import os
import awkward as ak
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tools.make_index_file import generate_html_indexes

colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#999999', '#ADFF2F', '#8B4513']

def make_plots(data_dict, output_dir):

    # Open the plot recipes`
    with open("tools/plot_recipes.json","r") as f:
        plot_recipes = json.load(f)

    # Grab the deault settings
    default = plot_recipes["default_settings"] 
 
    # Make the output directories
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "pLS"), exist_ok = True) 
    os.makedirs(os.path.join(output_dir, "LS"), exist_ok = True) 
    if default["plot_log_also"]:
        os.makedirs(os.path.join(output_dir, "pLS", "log"), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "LS", "log"), exist_ok=True)

    # Loop through each plot in the json (skipping default options)
    for plot_name, plot_cfg in plot_recipes.items():
        if plot_name == "default_settings":
            continue

        # Get the items that must be present in the plot
        title = plot_cfg["title"]
        objects = plot_cfg["objects"]
        legend = plot_cfg["legend"]

        # Get the rest of the options with default if needed
        xaxis_log = plot_cfg.get("xaxis_log", default["xaxis_log"])
        yaxis_log = plot_cfg.get("yaxis_log", False)
        plot_log_also = plot_cfg.get("plot_log_also", default["plot_log_also"])
        nbins = plot_cfg.get("nbins", default["nbins"])
        has_ratio = plot_cfg.get("has_ratio", default["has_ratio"])
        x_label_size = plot_cfg.get("x_label_size", default["x_label_size"])
        title_label_size = plot_cfg.get("title_label_size", default["title_label_size"])
        y_axis_label = plot_cfg.get("y_axis_label", default["y_axis_label"])
        # Special options for range
        if "range" in plot_cfg:
            plot_range = plot_cfg["range"]
        elif plot_name.endswith("_pt"):
            plot_range = default["pt_range"]
        elif plot_name.endswith("_phi"):
            plot_range = default["phi_range"]
        elif plot_name.endswith("_eta"):
            plot_range = default["eta_range"]
        else:
            raise Exception(f"Could not determine a valid range for plot: {plot_name}")
        # Special options for xlabel
        if "x_axis_label" in plot_cfg:
            xlabel = plot_cfg["x_axis_label"]
        elif plot_name.endswith("_pt"):
            xlabel = default["pt_xlabel"]
        elif plot_name.endswith("_phi"):
            xlabel = default["phi_xlabel"]
        elif plot_name.endswith("_eta"):
            xlabel = default["eta_xlabel"]
        else:
            raise Exception(f"Could not determine a valid x axis label for plot: {plot_name}")

        # We want plot_log_also to be False if yaxis_log is already true
        if yaxis_log:
            plot_log_also = False

        # Get the needed data and flatten it
        hist_data = [ak.to_numpy(ak.flatten(data_dict[obj])) for obj in objects]

        # Define the bin edges
        if xaxis_log:
            bin_edges = np.logspace(np.log10(max(plot_range[0],0.1)), np.log10(plot_range[1]), nbins + 1)
        else:
            bin_edges = np.linspace(plot_range[0], plot_range[1], nbins + 1)

        # make figure
        if has_ratio:
            fig = plt.figure(figsize=(8, 6))
            gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
            ax_main = fig.add_subplot(gs[0])
            ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)
        else:
            fig, ax_main = plt.subplots(figsize=(8, 5))

        # stacked histogram
        ax_main.hist(
            hist_data,
            bins=bin_edges,
            stacked=True,
            label=legend,
            color=colors[:len(hist_data)],
            histtype='stepfilled',
            edgecolor='black',
        )

        ax_main.set_ylabel(y_axis_label)
        ax_main.set_title(title, fontsize = title_label_size)
        ax_main.legend()

        if xaxis_log:
            ax_main.set_xscale('log')
        if plot_name.endswith("_pt"):
            ax_main.axvline(x=0.8, color='red', linestyle='--', linewidth=1.5, label='x = 0.8')

        if yaxis_log:
            ax_main.set_yscale('log')
            plot_log_also = False

        if has_ratio:
            raise Exception("Not Quite Ready for that yet")
        else:
            ax_main.set_xlabel(xlabel, fontsize = x_label_size) 
 
        # Determine subdirectory
        if "_ls_" in plot_name.lower():
            subdir = "LS"
        elif "_pls_" in plot_name.lower():
            subdir = "pLS"
        else:
            subdir = ""

        # Create filename and path
        base_filename = os.path.join(output_dir, subdir, plot_name) if subdir else os.path.join(output_dir, plot_name)

        # Save normal version
        fig.savefig(f"{base_filename}.png")
        fig.savefig(f"{base_filename}.pdf")
        plt.close(fig)
        print(f"Plot {plot_name} has been created!")


        # Create a log version of the plot
        if plot_log_also:

            # Re-make the plot with log y-axis
            if has_ratio:
                fig = plt.figure(figsize=(8, 6))
                gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
                ax_main = fig.add_subplot(gs[0])
                ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)
            else:
                fig, ax_main = plt.subplots(figsize=(8, 5))

            ax_main.hist(
                hist_data,
                bins=bin_edges,
                stacked=True,
                label=legend,
                color=colors[:len(hist_data)],
                histtype='stepfilled',
                edgecolor='black',
            )
            ax_main.set_ylabel(y_axis_label)
            ax_main.set_title(title, fontsize=title_label_size)
            ax_main.legend()
            ax_main.set_xlabel(xlabel, fontsize=x_label_size)

            if xaxis_log:
                ax_main.set_xscale('log')
            if plot_name.endswith("_pt"):
                ax_main.axvline(x=0.8, color='red', linestyle='--', linewidth=1.5, label='x = 0.8')

            ax_main.set_yscale('log')

            if has_ratio:
                raise Exception("Not Quite Ready for that yet")
            
            # Decide subdir for log
            if "_ls_" in plot_name.lower():
                log_dir = os.path.join(output_dir, "LS", "log")
            elif "_pls_" in plot_name.lower():
                log_dir = os.path.join(output_dir, "pLS", "log")
            else:
                log_dir = output_dir

            log_path = os.path.join(log_dir, plot_name + "_log")
            fig.savefig(f"{log_path}.png")
            fig.savefig(f"{log_path}.pdf")
            plt.close(fig)
            print(f"Plot {plot_name}_log has been created!")

    generate_html_indexes(output_dir)
    print("Index files generated!")
