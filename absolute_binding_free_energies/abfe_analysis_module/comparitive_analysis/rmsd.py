# Plot RMSD of protein accross desired proportion of trajectory for all
# stages and lambda windows and all runs

import matplotlib.pyplot as plt
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSD
from .restrained_dof import get_mda_universe
from ..get_data.dir_paths import get_dir_paths
from ..get_data.dir_paths import get_run_name
from ..save_data import mkdir_if_required


def plot_rmsd(leg, run_no, stage, lam_val, ax, percent_traj_use, selection):
    """Plot RMSD on supplied axis for given mda selection for specified 
    lam_val, stage, and run number. Reference frame taken as first frame not
    discarded.

    Args:
        leg (str): bound or free
        run_no (int): run number 
        stage (str): restrain, discharge, vanish
        lam_val (str): lambda window
        ax (matplotlib axis): axis on which to plot 
        percent_traj_use (float): Percentage of end of trajectory over which to 
        calculate RMSD. Reference frame taken as first frame not discarded.
        selection (str): atom selection in mda selection language
    """
    mobile = get_mda_universe(leg, run_no, stage, lam_val)

    # Get index of first frame to use and use this to get the reference
    n_frames = len(mobile.trajectory)
    first_frame_idx = round(n_frames - ((percent_traj_use/100)*n_frames)) # Index of first frame to be used for analysis
    print(f"Calculating RMSD for {leg} {stage} lambda: {lam_val} run:{run_no}")
    print(f"Reference taken from frame number: {first_frame_idx +1}")

    R = RMSD(mobile, select=selection, groupselections=[selection], ref_frame=first_frame_idx) # use specified frame of mobile as ref
    R.run()
    rmsd_results = R.results.rmsd.T
    time = [x/1000 for x in rmsd_results[1]] #convert to ns
    rmsd = rmsd_results[3]

    ax.plot(time[first_frame_idx+1:], rmsd[first_frame_idx+1:], label=f"Run {run_no}")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(r"RMSD ($\AA$)")
    
    # if lam_folders[i] == lam_folders[-1]:
    #     ax.set_xlabel("time (ps)")
    # 
    # if leg == "bound" and stage == "discharge":
    #      ax.set_ylabel(r"RMSD ($\AA$)")


def plot_rmsds(leg, runs, percent_traj_dict, selection):
    """Plot rmsds against time for all supplied runs, for 
    all lambda windows in all stages of given leg.

    Args:
        leg (str): bound or free
        runs (list): List of run numbers (ints)
        percent_traj_dict (dict): Specifies percentage of trajectories
        to use for given stage (of form {"stage":percent,...})
        selection (str): atom selection in mda selection language
    """
    print("###############################################################################################")
    print(f"Plotting RMSDs for {leg} leg and selection {selection}")

    paths = get_dir_paths(runs, leg)
    run_names = list(paths.keys())
    stages = list(paths[run_names[0]].keys())
    if "vanish" in stages:
        stages = ["vanish"] + stages # split vanish into two - easiest to add another vanish
    fig, _ = plt.subplots(1,1,figsize=(20,52), dpi=200)
    subfigs = fig.subfigures(1,1*len(stages)) # subfigs to allow different no plots in each column


    vanish_count = 0 # split vanish into two stages and count number of vanish stages already plotted

    for i, stage in enumerate(stages):
        lam_vals = paths[run_names[0]][stage]["lam_vals"]

        if stage == "vanish": # split in half as many windows
            median_idx = round(len(lam_vals)/2)
            if vanish_count == 0:
                lam_vals =lam_vals[:median_idx]
                vanish_count +=1
            else:
                lam_vals = lam_vals[median_idx:]

        if len(stages) == 1:
            subfig = subfigs
        else:
            subfig = subfigs[i]
        subfig.suptitle(f'{leg} {stage}')
        subfig_axs = subfig.subplots(len(lam_vals),1, sharex=True)
        for j, lam_val in enumerate(lam_vals):
            ax = subfig_axs[j]
            ax.title.set_text(f"$\lambda$ = {lam_val}")
            for k, run in enumerate(runs):
                plot_rmsd(leg,runs[k],stage,  # supply integer run no for sake of get_mda_universe
                lam_val,ax,percent_traj_dict[stage],selection)
            ax.legend(loc="best")
            # Keep only the bottom time 
            if j!= len(lam_vals) - 1:
                x_axis = ax.get_xaxis()
                x_axis.set_visible(False)
        
    mkdir_if_required("analysis/rmsd")
    fig.savefig(f'analysis/rmsd/{leg}_rmsd_{selection}.png')