# Based on convergence-plots.py (see https://github.com/michellab/MDM2-DG_paper)
# This produces convergence data for the stages specified (both overall and for
# individual stages). The resulting dictionary is of the form:
# {run:{stage:{cumtime1: {dg_tot:DG_tot, pmf:{lamval1: DG1, lamval2:DG2, ...}}, cumtime2:{...}, ...}}}

from itertools import starmap

from numpy import append
from . import dir_paths
import os, pickle
from multiprocessing import Pool, cpu_count
from ..save_data import mkdir_if_required


def truncate_simfile(in_file, out_file, start_time, end_time, 
                     nrg_freq, timestep, remove_heading=False):
    """Truncate a simfile between given start and end time.

    Args:
        in_file (str): Path to input file
        out_file (str): Path to input file 
        start_time (int): Start time in ns
        end_time (int): End time in ns
        nrg_freq (int): Number of steps between energy evaluations
        timestep (int): Timestep in ns
        remove_heading (bool): Whether or not to remove the heading
        before the data. Defaults to false.
    """
    with open(in_file, "r") as istream:
        with open(out_file, "w") as ostream:
            for line in istream.readlines():
                if line.startswith("#"):
                    if not remove_heading:
                        ostream.write(line)
                    continue
                elems = line.split()
                time = (float(elems[0])+nrg_freq)*timestep
                if (time < start_time):
                    continue
                if (time > end_time):
                    break
                ostream.write(line)


def do_mbar(lam_vals, input_dir, output_dir, start_time, end_time, nrg_freq, timestep):
    """Perform MBAR analysis for supplied lambda values for data between
    given start and end times. 

    Args:
        lam_vals (list): List of lambda values for which MBAR will be performed
        input_dir (str): Output directory for desired stage
        output_dir (str): Directory where output will be saved
        start_time (float): Start time in ns
        end_time (float): End time in ns
        nrg_freq (int): Number of steps between energy evaluations
        timestep (int): Timestep in ns

    Returns:
        (float): cumulative sampling time 
    """
    cumtime = 0
    delta_t = end_time - start_time

    for lam_val in lam_vals:
        os.system(f"mkdir {output_dir}/lambda-{lam_val}")
        truncate_simfile(f"{input_dir}/lambda-{lam_val}/simfile.dat",
                         f"{output_dir}/lambda-{lam_val}/simfile.dat",
                          start_time, end_time,
                          nrg_freq, timestep)
        cumtime += delta_t

    cmd = f"analyse_freenrg mbar \
         -i {output_dir}/lambda*/simfile.dat -p 100 --temperature 298.0 > {output_dir}/mbar.dat"
    os.system(cmd)
    print(f"MBAR analysis complete for cumulative sampling time per window {delta_t:.2f} ns")

    return cumtime


def read_mbar_data(mbar_file, lam_vals):
    """Read the total MBAR free energy estimate and MBAR
    PMF values

    Args:
        mbar_file (str): Path to .dat file to read
        lam_vals (list): List of lam vals in string format

    Returns:
        dg_tot (float): Total MBAR free energy for stage
        dg_sd (float): Standard deviation estimate from MBAR
        pmf (list): List of relative DG values from MBAR for each lambda value
        overlap (list): Matrix of overlap values (as a list of lists)
    """
    pmf = [] # List
    overlap = [] # List of lists (overlap matrix)
    no_lam_vals = len(lam_vals)

    with open(mbar_file,"r") as istream:
        lines = istream.readlines()

        for i,l in enumerate(lines):
            if l.startswith("#Overlap matrix"):
                for j in range(i+1, i+no_lam_vals+1):
                    overlap.append([float(x) for x in lines[j].strip().split(" ")])
            elif l.startswith("#PMF from MBAR in kcal/mol"):
                for j in range(i+1, i+no_lam_vals+1):
                    pmf.append(float(lines[j].strip().split(" ")[1]))
            elif l.startswith("#MBAR free energy difference"):
                dg_tot = float(lines[i+1].split()[0].strip(","))
                dg_conf_int = float(lines[i+1].split()[1].strip(","))

  #  for i,line in enumerate(overlap): # read in as str - convert to float
  #      for j, num in enumerate(line):
  #          overlap[i][j] = float(num)

    return dg_tot, dg_conf_int, pmf, overlap

def get_convergence_dict(leg="bound", run_nos=[1, 2, 3, 4, 5],
                        chunksize=0.05, nrg_freq=100, timestep=0.000004,
                        simtime = {"restrain": {"wind_len": 6, "discard": 1}, "discharge": {
                            "wind_len": 6, "discard": 1}, "vanish": {"wind_len": 8, "discard": 3},
                            "release": {"wind_len": 6, "discard": 1}, "unrigidify_lig": {
                            "wind_len": 6, "discard": 1},"unrigidify_prot": {
                            "wind_len": 6, "discard": 1}, "rigidify": {
                            "wind_len": 6, "discard": 1}, "release_2": {"wind_len": 2, "discard": 1}}
                         ):
    """Create dictionary of the form {run:{stage:{cumtime1: {dg_tot:DG_tot, 
    pmf:{lamval1: DG1, lamval2:DG2, ...}}, cumtime2:{...}, ...}}} to store all
    convergence data.

    Args:
        leg (string) : bound or free
        run_nos (list): Numbers of runs to include in analysis
        chunksize (float): ns between each calculation of energies
        nrg_freq (int): Number of steps between energy evaluations
        timestep (int): Timestep in ns
        simtime (dict): Lengths of simulations and lengths of inital
                        periods to discard as equilibration, in ns. 

    Returns:
        dict: convergence data
    """
    print("###############################################################################################")
    print("Calculating convergence data and obtaining convergence dictionary")

    paths = dir_paths.get_dir_paths(run_nos, leg)
    conv_dict = {}
    
    for run in paths.keys():
        conv_dict[run]={}
        for stage in paths[run].keys():
            conv_dict[run][stage]={}
            start_time = simtime[stage]["discard"]
            final_end_time = simtime[stage]["wind_len"]
            end_times = [(x*chunksize)+start_time for x in range(1, int((final_end_time-start_time)/chunksize)+1)]
            win_times = [x-start_time for x in end_times] # Cumulative time for single lam window
            lam_vals = paths[run][stage]["lam_vals"]
            input_dir = paths[run][stage]["output"] # output of simulations is input to do_mbar()

            print("###############################################################################################")
            print(f"MBAR analysis in progress for {run}, {stage}. Running in parallel over {cpu_count()} CPUs.")
            print("###############################################################################################")

            os.system("mkdir tmp")
            # create a list of tuples of the arguments to supply to starmap
            do_mbar_args = []
            for i, win_time in enumerate(win_times): # There will be a corresponding cumulative time, returned by do_mbar
                os.system(f"mkdir tmp/{win_time}")
                do_mbar_args.append((lam_vals, input_dir, f"./tmp/{win_time}",
                                      start_time, end_times[i], nrg_freq, timestep))

            # Carry out mbar analyses in parallel
            with Pool() as pool:
                cumtimes = pool.starmap(do_mbar, do_mbar_args)

            #extract data and delete dirs
            for i, win_time in enumerate(win_times): # There will be a corresponding cumulative time, returned by do_mbar
                cumtime = cumtimes[i]
                conv_dict[run][stage][cumtime] = {}
                dg_tot, _, pmf, _ = read_mbar_data(f"./tmp/{win_time}/mbar.dat",lam_vals) # throw away sd and overlap
                pmf_dict = {}
                for i, lam_val in enumerate(lam_vals):
                    pmf_dict[lam_val] = pmf[i]
                conv_dict[run][stage][cumtime]["dg_tot"] = dg_tot
                conv_dict[run][stage][cumtime]["pmf"] = pmf_dict
                
            os.system("rm -rf tmp")

    mkdir_if_required("analysis")
    dumpme = open("analysis/convergence_data.pickle","wb") # save data to file
    pickle.dump(conv_dict, dumpme)
    dumpme.close()

    return conv_dict

if __name__ == "__main__":
    get_convergence_dict()