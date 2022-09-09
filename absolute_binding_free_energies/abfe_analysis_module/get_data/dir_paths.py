# Takes input of list of run numbers. Returns dictionary containing all
# run directories, subdirectories, and a list of lam vals. Assumes directory
# structure base/bound/runXXX/stages and must be run from base.

import os
import re


def get_run_name(run_no, leg):
    """Finds top directory name for given run number and leg,
    ignoring comments on the end of the directory name.

    Args:
        run_no (int): Run number to find directory name for
        leg (str): Bound or free

    Returns:
        str: Directory name for specified run
    """

    run_no = str(run_no).zfill(3)  # format e.g. "001"
    dir_names = os.listdir(leg)
    # pattern match to ignore comments on dir names
    r = re.compile(f"run{run_no}")
    run_name = list(filter(r.match, dir_names))[0]

    return run_name


def get_file_name(path, pattern):
    """Finds file name containing pattern.

    Args:
        path (str): Path to desired folder
        pattern (str): Pattern in file name to match

    Returns:
        name: name of file matching pattern
    """

    file_names = os.listdir(path)
    # pattern match to ignore comments on names
    r = re.compile(pattern)
    try:
        name = list(filter(r.search, file_names))[0]
    except:
        print(f'Error: file matching "{pattern}" does not exist at "{path}"')
        name="no_file"

    return name


def get_subdir_names(directory):
    """Return list of subdirectory names.

    Args:
        directory (str): Path of directory

    Returns:
        list: Names of subdirectories
    """
    return [x for x in os.listdir(directory) if os.path.isdir(f"{directory}/{x}")]


def get_lam_vals(output_dir):
    """Get list of lambda directories and values based on directory 
    names in output directory.

    Args:
        output_dir (str): Path of output directory

    Returns:
        lam_dirs (list): Lambda directories
        lam_vals (list): Lambda values in string format
    """

    lam_dirs = [x for x in os.listdir(output_dir) if (("lambda" in x) and ("_" not in x))] #remove commented dirs
    lam_dirs.sort()
    lam_vals = [x[-5:] for x in lam_dirs]  # assume format "lambda-0.500"

    return lam_dirs, lam_vals


def get_dir_paths(run_nos, leg):
    """Return dictionary containing paths to input, output, and
    all lambda values. Also return list of lambda values.

    Args:
        run_nos (list): List of ints of run numbers to be included
        leg (str): Bound or free

    Returns:
        dict: Specifies paths
    """

    dir_paths = {}
    for run_no in run_nos:
        run_name = get_run_name(run_no, leg)
        dir_paths[run_name] = {}
        for stage in get_subdir_names(f"{leg}/{run_name}"):
            dir_paths[run_name][stage] = {}
            dir_paths[run_name][stage]["input"] = f"{leg}/{run_name}/{stage}/input"
            dir_paths[run_name][stage]["output"] = f"{leg}/{run_name}/{stage}/output"

            mbar_data_name = get_file_name(
                f"{leg}/{run_name}/{stage}/output", "MBAR")
            dir_paths[run_name][stage]["mbar_data"] = f"{leg}/{run_name}/{stage}/output/{mbar_data_name}"

            lam_dirs, lam_vals = get_lam_vals(
                f"{leg}/{run_name}/{stage}/output")
            dir_paths[run_name][stage]["lam_vals"] = lam_vals
            dir_paths[run_name][stage]["lam_paths"] = []
            for lam_dir in lam_dirs:
                dir_paths[run_name][stage]["lam_paths"].append(
                    f"{leg}/{run_name}/{stage}/output/{lam_dir}")

    return dir_paths


if __name__ == "__main__":
    import sys
    first_run = int(sys.argv[1])
    last_run = int(sys.argv[2])
    run_nos = [x for x in range(first_run, last_run+1)]
    leg = sys.argv[3]
    print(get_dir_paths(run_nos, leg))
