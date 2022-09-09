# Run complete set of default analyses 
import os

from .get_data import dir_paths
from .get_data import convergence_data
from .get_data import get_results
from .get_data import check_success
from .comparitive_analysis import compare_conv
from .comparitive_analysis import compare_pmfs
from .comparitive_analysis import overlap
from .comparitive_analysis import indiv_pmf_conv

def run_analysis(leg = "bound", run_nos=[1,2,3,4,5], restraint_type="Boresch", timestep=4, nrg_freq=100,
                 percent_traj_dict = {"restrain":83.33333333, "discharge":83.33333333, "vanish":62.5, "rigidify":83.33333333,
                     "unrigidify_lig":83.33333333, "unrigidify_prot":83.33333333, "release":83.333333},
                simtime = {"restrain": {"wind_len": 6, "discard": 1}, "discharge": {"wind_len": 6, "discard": 1}, "vanish": {"wind_len": 6, "discard": 1},
                          "release": {"wind_len": 6, "discard": 1}, "unrigidify_lig": {
                          "wind_len": 6, "discard": 1},"unrigidify_prot": {
                          "wind_len": 6, "discard": 1}, "rigidify": {
                          "wind_len": 6, "discard": 1}, "release_2": {"wind_len": 2, "discard": 1}}):
                
                    # For free leg, remember to change to:
                    #percent_traj_dict = {"discharge":83.33333333, "vanish":83.3333333}):
                    #simtime = {"discharge": {"wind_len": 6, "discard": 1}, "vanish": {"wind_len": 6, "discard": 1}} # ns

    """Run analysis of bound leg. Note that some global variables must be changed in convergence_data.py (this will
    be fixed in future). If the convergence data has already been generated, the analysis will start from there.

    Args:
        leg (str, optional): Does not currently allow 'free' as an option. Defaults to "bound".
        run_nos (list, optional): _description_. Defaults to [1,2,3,4,5].
        restraint_type (str, optional): Boresch, multiple_dist, or Cart. Defaults to "Boresch".
        timestep (int, optional): In fs. Defaults to 4.
        nrg_freq (int, optional): Steps between energy evaluations. Defaults to 100.
        percent_traj_dict (dict, optional): Percentage of trajectory to use for analysis for each stage.
        Defaults to {"restrain":83.33333333, "discharge":83.33333333, "vanish":62.5}.
        simtime (dict): Lengths of simulations and lengths of inital periods to discard as equilibration, in ns
                        , for generation of the convergence data. 
    """

    print("###############################################################################################")
    print(f"Analysing the {leg} leg for runs: {run_nos} and calculation type = {restraint_type}")
    print("Ensure you are in the base directory and the development version of biosimspace is activated")
    
    if leg == "free": # Different equilibration period and total sim time for vanish stage
        percent_traj_dict["vanish"] = 83.3333333
        simtime["vanish"] = {"wind_len": 6, "discard": 1}

    # Get stages
    paths = dir_paths.get_dir_paths(run_nos, leg)
    run_names = list(paths.keys())
    stages = list(paths[run_names[0]].keys())
    # Remove irrelevant entries from percent_traj_dict
    percent_traj_dict = {k:v for k,v in percent_traj_dict.items() if k in stages}
    simtime = {k:v for k,v in simtime.items() if k in stages}

    # Only calculate convergence data if this has not been done already
    if not os.path.isfile("analysis/convergence_data.pickle"):
        # Check if all simulations ran successfully
        if not check_success.check_success(run_nos=run_nos, leg=leg, verbose=True):
            raise Exception("Error: Not all windows completed successfully")
        convergence_data.get_convergence_dict(leg=leg, run_nos=run_nos, nrg_freq=nrg_freq,
                                              timestep=timestep/1000000, # Convert to ns
                                               simtime=simtime)

    if leg == "bound":
        get_results.write_results(leg, run_nos, restraint_type)
        compare_conv.plot_stages_conv("analysis/convergence_data.pickle", leg)
        compare_conv.plot_overall_conv("analysis/convergence_data.pickle", leg)
        compare_pmfs.plot_all_pmfs(run_nos, leg)
        overlap.plot_overlap_mats(leg, run_nos)
        indiv_pmf_conv.plot_pmfs_conv(leg, run_nos)

    elif leg == "free":
        get_results.write_results(leg, run_nos)
        compare_conv.plot_stages_conv("analysis/convergence_data.pickle", leg)
        compare_conv.plot_overall_conv("analysis/convergence_data.pickle", leg)
        compare_pmfs.plot_all_pmfs(run_nos, leg)
        overlap.plot_overlap_mats(leg, run_nos)
        indiv_pmf_conv.plot_pmfs_conv(leg, run_nos)

    print("###############################################################################################")
    print(f"Analysis of {leg} leg successfully completed!")