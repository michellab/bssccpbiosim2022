# Extracts all contributions to the free energy of binding for given leg (from MBAR)
# along with associated error calculated by propagation of the MBAR errors 
# (a large underestimation of the uncertainty calculated from multiple runs).

from cProfile import run
import numpy as np
import statsmodels.stats.api as sms

from .dir_paths import get_dir_paths
from ..save_data import mkdir_if_required
from .convergence_data import read_mbar_data
import scipy.stats as st


def get_lj_corr(input_file):
    """Finds LJ correction from .dat file

    Args:
        input_file (str): Input file

    Returns:
        : tuple of (LJ correction, standard deviation)
    """
    try:
        with open(input_file,'r') as file:
            lines = file.readlines()
        correction = float(lines[0].split()[2])
        conf_int= float(lines[0].split()[4])
        if (not np.isfinite(correction)) or (not np.isfinite(conf_int)):
            raise Exception()

    except:
        print(f"ERROR: Unable to read {input_file} or value in not finite. LJ correction has likely failed")
        correction = 0
        conf_int = 0
        
    return correction,conf_int


def get_restraint_cor(input_file, restraint_type):
    """Finds correction for releasing 
    restraints and returns correction as string

    Args:
        input_file (str): Input file
        restraint_type (str): Boresch, Cart, or multiple_dist

    Returns:
        energy: correction for releasing restraints
    """
    energy = 0
    conf_int = 0 # Include this for consistency with other functions, which return non-zero SD

    with open(input_file,'r') as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if restraint_type == "Boresch":
            if "correction for releasing Boresch restraints =" in line:
                energy=float(lines[i].split()[-3])
        if restraint_type == "Cart":
            if "correction for releasing Cartesian restraints =" in line:
                energy=float(lines[i].split()[-3])
        elif restraint_type == "multiple_dist":
            if "Free energy change upon removing the restraint" in line:
                energy=float(lines[i].split()[-2])

    if energy == 0:
        print(f"ERROR: Unable to read {input_file}. Correction has likely failed")
        energy = 0

    return energy, conf_int


def get_results(leg = "bound", run_nos = [1,2,3,4,5], restraint_type="Boresch"):
    """Get the MBAR free energy estimates and standard deviations
    associated with all stages in a given leg for the supplied runs

    Args:
        leg (str, optional): Free or bound. Defaults to "bound".
        run_nos (list, optional): List of int run numbers. Defaults to [1,2,3,4,5].
        restraint_type (str): Boresch, Cart, or multiple_dist

    Returns:
        dict: Dictionary of form {"run001":{"vanish":dg, "restrain":dg, "lj_corr":dg,...} ,...}
    """

    paths = get_dir_paths(run_nos, leg)
    results = {}
    
    for run_name in paths.keys():
        results[run_name]={}
        for stage in paths[run_name].keys():
            mbar_file = paths[run_name][stage]["mbar_data"]
            lam_vals = paths[run_name][stage]["lam_vals"]
            try:
                dg, conf_int, _, _ = read_mbar_data(mbar_file,lam_vals) # throw away PMF and overlap
            except Exception as e:
                print(f"Error: failure to write results for {run_name}, {stage}")
                raise e
            
            if stage in ["release", "unrigidify_lig", "unrigidify_prot"]: # Reverse sign of contribution
                dg *= -1
            results[run_name][stage] = (dg, conf_int)

            if stage == "vanish":
                output_dir = paths[run_name][stage]["output"]
                dg, conf_int = get_lj_corr(f"{output_dir}/freenrg-LJCOR.dat")
                results[run_name]["lj_corr"] = (dg, conf_int)
                if leg == "bound":
                    if restraint_type == "Boresch":
                        dg_ana, conf_int_ana = get_restraint_cor(f"{output_dir}/boresch_analytical_correction.dat",
                                                                restraint_type="Boresch")
                        results[run_name]["boresch_ana_corr"] = (dg_ana, conf_int_ana)
                        dg_semi, conf_int_semi = get_restraint_cor(f"{output_dir}/boresch_semi_ana_correction.dat",
                                                                    restraint_type="Boresch")
                        results[run_name]["boresch_semi-ana_corr"] = (dg_semi, conf_int_semi)
                    elif restraint_type == "multiple_dist":
                        dg_cor, conf_int_cor = get_restraint_cor(f"{output_dir}/standard-state-s-1-b-4-d-0.25-o-6.dat",
                                                                restraint_type="multiple_dist")
                        results[run_name]["multiple_dist_corr"] = (dg_cor, conf_int_cor)
                    elif restraint_type == "Cart":
                        dg_cor, conf_int_cor = get_restraint_cor(f"{output_dir}/cartesian_correction.dat",
                                                                restraint_type="Cart")
                        results[run_name]["cartesian_corr"] = (dg_cor, conf_int_cor)

                    # Symmetry corrections assume 298 K (RT = 0.592187)
                    results[run_name]["symm_corr_binding_sites_298"] = (0.65, 0) # Three-fold symmetry of binding sites (so RTln3)
                    #results[run_name]["symm_corr_phenol_298"] = (0.41, 0) # Rotation of phenol hindered in binding site (so RTln2)
            
        
        dg_tot = sum([val[0] for key, val in results[run_name].items() if key != "boresch_ana_corr"]) # Energy is first value in tuple.
                                                                                                      # Avoid adding two Boresch corrections.
        var_tot = sum([x[1]**2 for x in results[run_name].values()])
        ci_tot = np.sqrt(var_tot)
        results[run_name]["dg_tot"] = (dg_tot, ci_tot) 
                    
    return results


def write_results_indiv(results):
    """Write results for individual runs to text file.

    Args:
        results (dict): Dictionary of results (possibly for several runs)
    """
    for run_name in results.keys():
        mkdir_if_required("analysis")
        mkdir_if_required("analysis/results")
        mkdir_if_required(f"analysis/results/{run_name}")

        with open(f"analysis/results/{run_name}/{run_name}_results.txt", "wt") as f:
            f.write("Free energy estimates from MBAR:\n")
            f.write("(Uncertainties are 95 % C.I.s derived from MBAR for single run)\n\n")
            for contribution in results[run_name].keys():
                dg = results[run_name][contribution][0]
                sd = results[run_name][contribution][1]
                f.write(f"{contribution}: {dg:.2f} +/- {sd:.2f} kcal/mol\n")


def write_results_overall(results):
    """Write summary results to text file. Calculate
    95 % C.I. for all stages.

    Args:
        results (dict): Results dictionary
    """
    mkdir_if_required("analysis")
    mkdir_if_required("analysis/results")

    tot_dict = {}
    conf_ints = [] # Store C.I.s for all contributions and use these to work out C.I. for dg_tot
    run_names = list(results.keys())
    for contribution in results[run_names[0]]: # Use first run to check contributions to free energy
        tot_dict[contribution]={}
        tot_dict[contribution]["values"]=[]
        for run_name in run_names:
            tot_dict[contribution]["values"].append(results[run_name][contribution][0]) # Ignore SD from individual runs

        tot_dict[contribution]["values"] =np.array([x for x in tot_dict[contribution]["values"] if np.isfinite(x)])
        vals = tot_dict[contribution]["values"]
        tot_dict[contribution]["dg"] = vals.mean()

        if contribution != "dg_tot":
            conf_int = st.t.interval(0.95, len(vals)-1, loc=np.mean(vals), scale=st.sem(vals))
            positive_conf_int = vals.mean() - conf_int[0] # Because C.I. returned as tuple (min, max)
            tot_dict[contribution]["95% C.I."] = positive_conf_int
            conf_ints.append(positive_conf_int)

        elif contribution == "dg_tot":
        # Get C.I. by adding C.I.s in quadrature for individual contributions, rather than using C.I.s based on dg_tots
        # from all runs, as this results in the loss of information. NOTE: This relies on dg_tot being the last contribution
            conf_ints = np.array(conf_ints)
            filtered_conf_ints = np.array([x for x in conf_ints if np.isfinite(x)])
            tot_dict[contribution]["95% C.I."] = np.sqrt(sum(filtered_conf_ints**2))


    with open(f"analysis/results/results.txt", "wt") as f:
        f.write("Overall free energy estimates from MBAR:\n")
        f.write(f"(Uncertainty estimates are 95 % C.I.s derived from differences between {len(run_names)} replicate runs)\n\n")
        for contribution in tot_dict.keys():
            dg = tot_dict[contribution]["dg"]
            conf = tot_dict[contribution]["95% C.I."]
            f.write(f"{contribution}: {dg:.2f} +/- {conf:.2f} kcal/mol\n")


def write_results(leg="bound", run_nos = [1,2,3,4,5], restraint_type="Boresch"):
    """Retrieve results for given leg for all supplied runs.
    Save individual summary and overall summary with 95 % C.I.s

    Args:
        leg (str, optional): Bound or free. Defaults to "bound".
        run_nos (list, optional): List of run numbers (ints). Defaults to [1,2,3,4,5].
        restraint_type (str): Boresch, Cart, or multiple_dist
    """
    print("###############################################################################################")
    print("Writing results")
    results = get_results(leg, run_nos, restraint_type)
    write_results_indiv(results)
    write_results_overall(results)