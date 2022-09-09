"""Check that all calculations completed successfully; list failed windows and reason if not"""

from .dir_paths import get_dir_paths
import glob

def check_success(run_nos = [1,2,3,4,5], leg="bound", verbose=True):
    """Check that all lambda windows completed successfully; list
    failures if not.

    Args:
        run_nos (list, optional): Run numbers to check. Defaults to [1,2,3,4,5].
        leg (str, optional): bound or free. Defaults to bound.
    """
    print("###############################################################################################")
    print("Checking for successful completion of simulations")
    failed_windows = []
    no_instability_failures = 0
    dir_paths = get_dir_paths(run_nos, leg)
    for run in dir_paths:
        for leg in dir_paths[run]:
            output_dir = dir_paths[run][leg]["output"]
            slurm_logs = glob.glob(f"{output_dir}/somd-array-gpu*", recursive=False)
            for log in slurm_logs:
                lam = "Undefined"
                with open(log, "rt") as f:
                    success = False
                    lines = f.readlines()
                    no_lines = len(lines)
                    for l in lines:
                        if l.startswith("lambda is"):
                            lam = l.split()[-1]
                            break
                    for l in lines[no_lines - 10:]: # To deal with massive debug == True outputs
                        if l.startswith("lambda is"):
                            lam = l.split()[-1]
                        if "NaN or Inf has been generated along the simulation" in l:
                            no_instability_failures += 1
                        if "RuntimeError: Particle coordinate is nan" in l:
                            no_instability_failures += 1
                        if l.startswith("Simulation took"):
                            success = True
                            break
                    if not success:
                        failure_reason = (l for l in [lines[no_lines - 2], lines[no_lines - 1]])
                        failure_reason = ', '.join(failure_reason)
                        failed_windows.append((f"{run} {leg} {lam}", failure_reason))

    if failed_windows == []:
        print("All windows ran successfully")
        return True
    else:
        print("The following windows failed:")
        if verbose:
            for fail in sorted(failed_windows):
                print(f"{fail[0]}: {fail[1]}")
        else:
            for fail in sorted(failed_windows):
                print(f"{fail[0]}")
        print(f"No failed windows = {len(failed_windows)}")
        print(f"Failures due to instabilities = {no_instability_failures}")

        return False

