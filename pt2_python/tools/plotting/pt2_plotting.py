import os
import awkward as ak
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tools.plotting.make_index_file import generate_html_indexes

colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#999999', '#ADFF2F', '#8B4513']

def make_plots(data_dict, output_dir):

    # Load plot recipes json
    with open("tools/plotting/plot_recipes.json", "r") as f:
        plot_recipes = json.load(f)

    # Get the default settings
    default = plot_recipes["default_settings"]

    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "pLS"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "LS"), exist_ok=True)
    if default["plot_log_also"]:
        os.makedirs(os.path.join(output_dir, "pLS", "log"), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "LS", "log"), exist_ok=True)

    for plot_name, plot_cfg in plot_recipes.items():
        if plot_name == "default_settings":
            continue
        if not "lst_objects" in plot_cfg:
            raise Exception(f"{plot_name} has no defined objects!")
        for lst_obj in plot_cfg["lst_objects"]:
            if not "values" in plot_cfg:
                raise Exception(f"{plot_name} has no values!")
            for val in plot_cfg["values"]:
                if plot_name.startswith("simIdx"):
                    ref_obj = "ls" if lst_obj.startswith("ls") else "pls"
                    other_obj = "pls" if ref_obj == "ls" else "ls"
                    used_unused = "Used" if "used" in lst_obj else "Unused"
                    title = f"{used_unused} {ref_obj}  " + plot_cfg["title"]
                    legend = [f"Real {ref_obj}: {leg} {other_obj}" for leg in plot_cfg["legend"][:-1]]
                    legend.append(f"{plot_cfg['legend'][-1]} {ref_obj}")
                    objects = [f"{lst_obj}_{base}_{other_obj}_{val}" for base in plot_cfg["objects_base"][:-1]]
                    objects.append(f"{lst_obj}_{plot_cfg['objects_base'][-1]}_{val}")

                else:
                    objects = [f"{lst_obj}_{base}_{val}" for base in plot_cfg["objects_base"]]
                    title = plot_cfg["title"] + f" {lst_obj} {val}"
                    legend = [f"{leg} {lst_obj}" for leg in plot_cfg["legend"]]

                xaxis_log = True if val == "pt" else plot_cfg.get("xaxis_log", default["xaxis_log"])
                yaxis_log = plot_cfg.get("yaxis_log", False)
                plot_log_also = plot_cfg.get("plot_log_also", default["plot_log_also"])
                nbins = plot_cfg.get("nbins", default["nbins"])
                has_ratio = plot_cfg.get("has_ratio", default["has_ratio"])
                x_label_size = plot_cfg.get("x_label_size", default["x_label_size"])
                title_label_size = plot_cfg.get("title_label_size", default["title_label_size"])
                y_axis_label = plot_cfg.get("y_axis_label", default["y_axis_label"])

                # Determine axis range
                if "range" in plot_cfg:
                    plot_range = plot_cfg["range"]
                elif val == "pt":
                    plot_range = default["pt_range"]
                elif val == "phi":
                    plot_range = default["phi_range"]
                elif val == "eta":
                    plot_range = default["eta_range"]
                else:
                    raise Exception(f"Could not determine a valid range for plot: {plot_name}")

                # Determine axis label
                if "x_axis_label" in plot_cfg:
                    xlabel = plot_cfg["x_axis_label"]
                elif val == "pt":
                    xlabel = default["pt_xlabel"]
                elif val == "phi":
                    xlabel = default["phi_xlabel"]
                elif val == "eta":
                    xlabel = default["eta_xlabel"]
                else:
                    raise Exception(f"Could not determine a valid x axis label for plot: {plot_name}")

                if yaxis_log:
                    plot_log_also = False

                hist_data = [ak.to_numpy(ak.flatten(data_dict[obj])) for obj in objects]

                if xaxis_log:
                    bin_edges = np.logspace(np.log10(max(plot_range[0], 0.1)), np.log10(plot_range[1]), nbins + 1)
                else:
                    bin_edges = np.linspace(plot_range[0], plot_range[1], nbins + 1)

                def draw_legend(ax):
                    if default["plot_legend_onside"]:
                        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=10)
                    else:
                        ax.legend()

                def get_output_path(is_log_version=False):
                    if lst_obj == "ls" or lst_obj.startswith("ls"):
                        subdir = "LS"
                    elif lst_obj == "pls" or lst_obj.startswith("pls"):
                        subdir = "pLS"
                    else:
                        subdir = ""

                    if is_log_version:
                        base_dir = os.path.join(output_dir, subdir, "log") if subdir else os.path.join(output_dir, "log")
                    else:
                        base_dir = os.path.join(output_dir, subdir) if subdir else output_dir

                    os.makedirs(base_dir, exist_ok=True)
                    if plot_name.startswith("simIdx"):
                        return os.path.join(base_dir, f"{lst_obj}_{plot_name}_{other_obj}_{val}" + ("_log" if is_log_version else ""))
                    else:
                        return os.path.join(base_dir, f"{lst_obj}_{plot_name}_{val}" + ("_log" if is_log_version else ""))

                def create_and_save_plot(is_log_version=False):
                    fig, ax_main = plt.subplots(figsize=(10, 5) if default["plot_legend_onside"] else (8, 5))

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
                    draw_legend(ax_main)
                    ax_main.set_xlabel(xlabel, fontsize=x_label_size)

                    if xaxis_log:
                        ax_main.set_xscale('log')
                    if val == "pt":
                        ax_main.axvline(x=0.8, color='red', linestyle='--', linewidth=1.5, label='x = 0.8')
                    if is_log_version or yaxis_log:
                        ax_main.set_yscale('log')

                    path = get_output_path(is_log_version)
                    fig.savefig(f"{path}.png", bbox_inches='tight')
                    fig.savefig(f"{path}.pdf", bbox_inches='tight')
                    plt.close(fig)
                    print(f"Plot {path} has been created!")

                if has_ratio:
                    raise Exception("Not Quite Ready for that yet")
                else:
                    create_and_save_plot(is_log_version=False)

                if plot_log_also:
                    create_and_save_plot(is_log_version=True)


    generate_html_indexes(output_dir)
    print("Index files generated!")
