"""Compare the results of two different sets of runs and test
significance of differences"""

from asyncio import run
from cmath import nan
from .get_results import get_results
import os
import scipy.stats as st
import numpy as np


def get_combined_results(calculation_1_path, calculation_2_path, leg = "bound", run_nos_calc_1 = [1,2,3,4,5],
                        run_nos_calc_2 = [1,2,3,4,5], restraint_type_calc_1="Boresch", restraint_type_calc_2="Boresch"):
    """Returns dictionary containing the results for each contribution to the free energy for each
    replicate for the calculations given by each of the supplied paths.

    Args:
        calculation_1_path (str): Path to base directory of first set of calculations to be compared
        calculation_2_path (str): Path to base directory of second set of calculations to be compared
        leg (str, optional): Bound or free. Defaults to "bound".
        run_nos_calc_1 (list, optional): Replicate runs to be included for run 1. Defaults to [1,2,3,4,5].
        run_nos_calc_2 (list, optional): Replicate runs to be included for run 2. Defaults to [1,2,3,4,5].
        restraint_type_calc_1(str, optional): Boresch, multiple_dist, or Cart. Defaults to Boresch.
        restraint_type_calc_2(str, optional): Boresch, multiple_dist, or Cart. Defaults to Boresch.

    Returns:
        dict: combined results
    """
    cwd = os.getcwd() # Get current working dir so we can return
    combined_results = {}

    for i, path_type in enumerate(zip([calculation_1_path, calculation_2_path], [restraint_type_calc_1, restraint_type_calc_2])):
        os.chdir(path_type[0])
        dir_name = os.path.split(os.getcwd())[1] # Use dir name to name the dictionary - "." uninformative, for example
        if i == 0:
            results_dict = get_results(leg, run_nos_calc_1, restraint_type=path_type[1])
        else:
            results_dict = get_results(leg, run_nos_calc_2, restraint_type=path_type[1])
        run_list = list(results_dict.keys())
        overall_results_dict = {}
        for contribution in results_dict[run_list[0]]:
            if contribution[-4:] == "corr" and contribution != "lj_corr":
                overall_results_dict["restraint_corr"] = [] # Name all corrections the same
            else:
                overall_results_dict[contribution] =[]

        for run in results_dict:
            results_subdict = results_dict[run]
            for contribution in results_subdict:
                if contribution[-4:] == "corr" and contribution != "lj_corr": # Make sure corrections have the same name so code below will work
                    overall_results_dict["restraint_corr"].append(results_subdict[contribution][0]) # As dg first in tuple
                else:
                    overall_results_dict[contribution].append(results_subdict[contribution][0]) # As dg first in tuple

        combined_results[dir_name] = overall_results_dict
        os.chdir(cwd)

    # Ensure all contributions are present in both sub-dictionaries
    dict1_id, dict2_id = combined_results.keys()
    for contribution in combined_results[dict1_id]:
        if contribution not in combined_results[dict2_id]:
            combined_results[dict2_id][contribution]=[]
    for contribution in combined_results[dict2_id]:
        if contribution not in combined_results[dict1_id]:
            combined_results[dict1_id][contribution]=[]

    return combined_results


def independent_ttest(mean1, sd1, n1, mean2, sd2, n2):
    """Perform independent t-test assuming equal variance.

    Args:
        mean1 (float): Mean of set 1
        sd1 (float): Standard deviation of set 1
        n1 (int): Number of samples in set 1
        mean2 (float): Mean of set 2
        sd2 (float): Standard deviation of set 2
        n2 (int): Number of samples in set 2

    Returns:
        float: p value
    """
    df = n1 + n2 -2

    # calculate pooled standard deviation (assumes the same)
    numerator_pooled = (n1 - 1)*sd1**2 + (n2 - 1)*sd2**2
    pooled_sd = np.sqrt(numerator_pooled/df)

    # calculate t
    numerator_t = mean1 - mean2
    denominator_t = pooled_sd*(np.sqrt((1/n1)+(1/n2)))
    t = numerator_t/denominator_t

    # p value
    p = st.t.sf(np.abs(t), df)*2  # two-sided pvalue = Prob(abs(t)>tt)

    return p
    

def write_sig_diff(calculation_1_path, calculation_2_path, leg = "bound", run_nos_calc_1 = [1,2,3,4,5], 
                   run_nos_calc_2 = [1,2,3,4,5], offset_calc_2=0, restraint_type_calc_1="Boresch", restraint_type_calc_2="Boresch"):
    """Carry out unpaired t-test assuming equal variances and write out p values. Note that the variance
    for dg_tot is calculated from the variances of the contributions, rather than from variation in dg_tot
    between runs, to avoid losing information.

    Args:
        calculation_1_path (str): Path to base directory of first set of calculations to be compared
        calculation_2_path (str): Path to base directory of second set of calculations to be compared
        leg (str, optional): Bound or free. Defaults to "bound".
        run_nos_calc_1 (list, optional): Replicate runs to be included for run 1. Defaults to [1,2,3,4,5].
        run_nos_calc_2 (list, optional): Replicate runs to be included for run 2. Defaults to [1,2,3,4,5].
        offset_calc_2 (float): Add an offset to the total free energy for calculation 2 relative to
        calculation 1. This is needed in the case of different symmetry corrections between runs.
        restraint_type_calc_1 (str, optional): Boresch, multiple_dist, or Cart. Defaults to Boresch.
        restraint_type_calc_2 (str, optional): Boresch, multiple_dist, or Cart. Defaults to Boresch.
    """

    combined_results = get_combined_results(calculation_1_path, calculation_2_path, leg, run_nos_calc_1, 
                                            run_nos_calc_2, restraint_type_calc_1, restraint_type_calc_2)

    dir1, dir2 = combined_results
    with open(f"{dir1}_{dir2}_sig_diff.txt", "wt") as f:
        f.write(f"{dir1} compared to {dir2}, unpaired t-test assuming equal variances\n")
        f.write(f"Calculation 2 overall results offset by {offset_calc_2} kcal mol-1.\n\n")

        stds_1 = []
        stds_2 = []

        for contribution in combined_results[dir1]:
            results_1 = np.array([x for x in combined_results[dir1][contribution] if np.isfinite(x)])
            results_2 = np.array([x for x in combined_results[dir2][contribution] if np.isfinite(x)])

            if contribution != "dg_tot":
                stds_1.append(results_1.std(ddof =1)) # Otherwise divides by n rather than n -1
                stds_2.append(results_2.std(ddof =1)) # Otherwise divides by n rather than n -1
                p = st.ttest_ind(results_1, results_2)[1]

            if contribution == "dg_tot":
                # Because we get the C.I. based on the C.I.s of the components and not from the overall results 
                # from each run (to preserve information), we have to calculate p manually from the variances
                # NOTE: This assumes that dg_tot is the last contribution
                results_2 += offset_calc_2
                stds_1 = np.array([x for x in stds_1 if np.isfinite(x) and  x!= 0])
                stds_2 = np.array([x for x in stds_2 if np.isfinite(x) and x!= 0])
                n1 = len(run_nos_calc_1)
                n2 = len(run_nos_calc_2)
                #filtered_stds_1 = [x for x in stds1 if x != nan] 
                #filtered_stds_2 = [x for x in stds2 if x != nan] 
                sd_tot_1 = np.sqrt(sum(stds_1**2))
                sd_tot_2 = np.sqrt(sum(stds_2**2))
                p = independent_ttest(results_1.mean(), sd_tot_1, n1, results_2.mean(), sd_tot_2, n2)

            line = f"{contribution} = "
            line += f"{results_1.mean():.3f}, {results_2.mean():.3f} "

            line += f", t-test p = {p:.3f}"
            if p < 0.05:
                line += ", SIGNIFICANT DIFFERENCE AT 95 % CONFIDENCE\n"
            else:
                line += "\n"
            f.write(line)