# Plot convergence of the pmf with time for each individual run

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import pickle

from ..get_data.dir_paths import get_dir_paths
from ..save_data import mkdir_if_required

def plot_pmf_conv(ax, leg, stage, run_name, conv_dict):
    """Plot convergence of PMF with cumulative sampling time
    for given leg, stage, and run, on supplied axis.

    Args:
        ax (matplotlib axis): Axis on which to plot
        leg (str): bound or free
        stage (str): e.g. vanish
        run_name (str): e.g. run001
        conv_dict (dict): The unpickled convergence data

    Returns:
        mapper: Colour mapper for cumulative sampling time per window.
        Allows this to be placed in selected axes in figure, rather than
        all.

    """
    print(f"Plotting convergence of PMF for run {run_name}, {leg} {stage}")
    cumtime_dict = conv_dict[run_name][stage]
    cum_times = list(cumtime_dict.keys())
    lam_vals_str = list(cumtime_dict[cum_times[0]]["pmf"].keys())
    lam_vals = [float(x) for x in lam_vals_str]
    n_lam_vals = len(lam_vals)
    cum_times_per_wind = [x/n_lam_vals for x in cum_times] # Because cum_times are cumulative over all windows

    # create mapper to map sampling time to colour
    norm = colors.Normalize(vmin=cum_times[0], vmax=cum_times_per_wind[-1], clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.brg)

    for i, cum_time in enumerate(cum_times):
        pmf = list(cumtime_dict[cum_time]["pmf"].values())
        print(cum_times_per_wind[i])
        ax.plot(lam_vals, pmf, c=mapper.to_rgba(cum_times_per_wind[i]), alpha=0.2, lw=0.5)
    ax.set_title(f'{run_name} {leg} {stage}')
    ax.set_xlabel(r'$\lambda$')
    ax.set_ylabel('$\Delta \it{G}$ / kcal.mol$^-$$^1$')

    return mapper


def plot_pmfs_conv(leg="bound", run_nos=[1,2,3,4,5], pickled_data="analysis/convergence_data.pickle"):
    """Plot convergence of PMFs for individual stages of individual runs in grid.

    Args:
        leg (str, optional): Bound or free. Defaults to "bound".
        run_nos (list, optional): List of run numbers to plot. Defaults to [1,2,3,4,5].
        pickled_data (str, optional): Path to pickled data file. Defaults to "analysis/convergence_data.pickle".
    """
    print("###############################################################################################")
    print("Plotting convergence of all PMFs individually")

    with open(pickled_data, "rb") as istream:
        conv_dict = pickle.load(istream)

    paths = get_dir_paths(run_nos, leg)
    run_names = list(paths.keys())
    stages = list(paths[run_names[0]].keys())
    n_runs = len(run_names)
    n_stages = len(stages)

    fig, axs = plt.subplots(n_runs, n_stages, figsize=(4*n_stages, 4*n_runs),dpi=1000)

    for i, run_name in enumerate(run_names):
        for j, stage in enumerate(stages):
            if n_stages ==1:
                ax = axs[i]
            if n_runs ==1:
                ax = axs[j]
            else:
                ax = axs[i,j]
            mapper = plot_pmf_conv(ax, leg, stage, run_name, conv_dict)
            if j == n_stages-1:
                fig.colorbar(mapper, ax=ax).set_label('Cumulative Sampling Time per Window / ns')
        
    fig.tight_layout()
    mkdir_if_required("analysis/individual")
    fig.savefig(f"analysis/individual/{leg}_individual_pmf_convergence.png")