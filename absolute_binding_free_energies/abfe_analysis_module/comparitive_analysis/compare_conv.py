# Plot comparitive convergence (overall and for individual stages) for
# all runs in conv_dict, which is of the form {run:{stage:{cumtime1:
# {dg_tot:DG_tot, pmf:{lamval1: DG1, lamval2:DG2, ...}}, cumtime2:{...}, ...}}} 
# Assumes sample size of 5

import matplotlib.pyplot as plt
import pickle
import numpy as np
from ..save_data import mkdir_if_required
import scipy.stats as st


def plot_conv(ax, leg, stage, x_vals, y_vals_list, x_label, y_label):
    """Plot convergence with cumulative sampling time 
    for a given stage.

    Args:
        ax (ax): Axis on which to plot
        leg (str): Bound or free
        stage (stage): Restrain, discharge, or vanish
        x_vals (list): List of values to be plotted. Can be string as well as int or float format
        y_vals_list (list): List of lists values to be plotted (e.g.PMFS)
        x_label (str): x-axis label
        y_label (str): y-axis label
    """
    x_vals = np.array(x_vals).astype(float)
    y_vals = np.array(y_vals_list).astype(float) # No problem even if input is already array
    y_avg = np.mean(y_vals, axis=0)
    conf_int = st.t.interval(0.95, len(y_vals[:,0])-1, loc=y_avg, scale=st.sem(y_vals,axis=0)) # 95 % C.I.

    #plt.xticks(np.linspace(min(x_vals), max(x_vals),6))
    ax.plot(x_vals,y_avg)
    for i, entry in enumerate(y_vals):
        ax.plot(x_vals,entry, linestyle='dashed', label=f"run {i+1}")
    ax.set_xscale("linear")
    ax.set_title(f'{leg} {stage}')
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    for tick in ax.get_xticklabels():
        tick.set_rotation(-45)

    ax.legend()
    ax.fill_between(x_vals, conf_int[0], conf_int[1], alpha=0.5, facecolor='#ffa500')


def plot_stages_conv(pickled_data, leg):
    """Plot convergence accross all runs for individual stages.

    Args:
        pickled_data (str): Path to pickled convergence dictionary data
        leg (str): Bound or free
    """
    print("###############################################################################################")
    print("Plotting convergence for individual stages")

    with open(pickled_data, "rb") as istream:
        conv_dict = pickle.load(istream)

    runs = list(conv_dict.keys())
    stages = list(conv_dict[runs[0]].keys())
    no_stages = len(stages)

    fig, axs = plt.subplots(1, no_stages, figsize = (4*no_stages,4), dpi = 1000)
    for i, stage in enumerate(stages):
        cumtimes = list(conv_dict[runs[0]][stage].keys())
        fr_nrgs = []
        for run in runs:
            fr_nrg_run = []
            for cumtime in cumtimes:
                fr_nrg_run.append(conv_dict[run][stage][cumtime]["dg_tot"])
            fr_nrgs.append(fr_nrg_run)

        if len(stages) == 1:
            plot_conv(axs, leg, stage, cumtimes, fr_nrgs, 'Cumulative sampling time per run / ns', '$\Delta \it{G}$ / kcal.mol$^-$$^1$')
        else:
            plot_conv(axs[i], leg, stage, cumtimes, fr_nrgs, 'Cumulative sampling time per run / ns', '$\Delta \it{G}$ / kcal.mol$^-$$^1$')
        

    fig.tight_layout()
    mkdir_if_required("analysis/overall_convergence")
    fig.savefig(f"analysis/overall_convergence/{leg}_stages_joint_convergence.png")


def plot_overall_conv(pickled_data, leg, exclude=[]):
    """Plot overall convergence accross all runs for given leg.

    Args:
        pickled_data (str): Path to pickled convergence dictionary data
        leg (str): Bound or free
        exclude (list): List of stages to exclude from the overall convergence plot.
    """
    print("###############################################################################################")
    print("Plotting overall convergence")

    with open(pickled_data, "rb") as istream:
        conv_dict = pickle.load(istream)

    runs = list(conv_dict.keys())
    stages = list(conv_dict[runs[0]].keys())
    fig, ax = plt.subplots(1, 1, figsize = (4,4), dpi = 1000)

    tot_cumtimes = []
    tot_fr_nrgs = []

    for stage in stages:
        if stage in exclude:
            pass
        else:
            cumtimes = list(conv_dict[runs[0]][stage].keys())
            tot_cumtimes.append(cumtimes)
            fr_nrgs = []
            for run in runs:
                fr_nrg_run = []
                for cumtime in cumtimes:
                    if stage in ["release", "unrigidify_lig", "unrigidify_prot"]: # Reverse sign of contribution
                        fr_nrg_run.append(-conv_dict[run][stage][cumtime]["dg_tot"])
                    else:
                        fr_nrg_run.append(conv_dict[run][stage][cumtime]["dg_tot"])
                fr_nrgs.append(fr_nrg_run)
            tot_fr_nrgs.append(fr_nrgs)
        
    tot_cumtimes = np.sum(np.array(tot_cumtimes), axis=0)
    tot_fr_nrgs = np.sum(np.array(tot_fr_nrgs), axis=0)

    plot_conv(ax, leg, "overall", tot_cumtimes, tot_fr_nrgs, 'Cumulative sampling time per run / ns', '$\Delta \it{G}$ / kcal.mol$^-$$^1$')
    fig.tight_layout()
    mkdir_if_required("analysis/overall_convergence")
    if not exclude:
        fig_name = f"analysis/overall_convergence/{leg}_overall_convergence.png"
    else:
        excluded = "_".join(exclude)
        fig_name = f"analysis/overall_convergence/{leg}_overall_convergence_exclude_{excluded}.png"
    fig.savefig(fig_name)


if __name__ == "__main__":
    import sys
    data = sys.argv[1]
    plot_stages_conv(data)
    plot_overall_conv(data)