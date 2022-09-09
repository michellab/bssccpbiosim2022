# Plot all overlap matrices for a given run

import numpy as np
import matplotlib.pylab as plt
import seaborn as sbn

from ..get_data.dir_paths import get_dir_paths
from ..save_data import mkdir_if_required
from ..get_data.convergence_data import read_mbar_data


def plot_overlap_mat(ax, leg, stage, run_name, mbar_file, lam_vals):
    """Plot overlap
    for given leg, stage, and run, on supplied axis.

    Args:
        ax (matplotlib axis): Axis on which to plot
        leg (str): bound or free
        stage (str): e.g. vanish
        run_name (str): e.g. run001
        mbar_file (dict): Path to MBAR data file
    """
    print(f"Plotting overlap matrix for run {run_name}, {leg} {stage}")
    
    _,_,_,overlap = read_mbar_data(mbar_file, lam_vals) # Throw away dg, etc.
    overlap = np.array(overlap)
    sbn.heatmap(overlap, ax=ax, square=True, linewidths=.5).figure
    ax.set_title(f'{run_name} {leg} {stage}')


def plot_overlap_mats(leg="bound", run_nos=[1,2,3,4,5]):
    """Plot all overlap matrices for given leg and given runs.

    Args:
        leg (str, optional): Bound or free. Defaults to "bound".
        run_nos (list, optional): Run numbers for runs to plot. Defaults to [1,2,3,4,5].
    """
    print("###############################################################################################")
    print("Plotting all overlap matrices")
    
    paths = get_dir_paths(run_nos, leg)
    run_names = list(paths.keys())
    stages = list(paths[run_names[0]].keys())
    n_runs = len(run_names)
    n_stages = len(stages)

    fig, axs = plt.subplots(n_runs, n_stages, figsize=(4*n_stages, 4*n_runs),dpi=1000)

    for i, run_name in enumerate(run_names):
        for j, stage in enumerate(stages):
            mbar_file = paths[run_name][stage]["mbar_data"]
            lam_vals_str = paths[run_name][stage]["lam_vals"]
            lam_vals = [float(x) for x in lam_vals_str]
            if n_stages == 1:
                ax = axs[i]
            if n_runs == 1:
                ax = axs[j]
            else:
                ax = axs[i,j]
            plot_overlap_mat(ax, leg, stage, run_name, mbar_file, lam_vals)
        
    fig.tight_layout()
    mkdir_if_required("analysis/individual")
    fig.savefig(f"analysis/individual/{leg}_overlap.png")