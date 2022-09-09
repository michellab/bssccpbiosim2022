# Functions to get histograms of Boresch degrees of freedom over
# trajectories for specified value of lambda

from math import ceil
import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from ..get_data import dir_paths
from ast import literal_eval
from ..save_data import mkdir_if_required


def get_mda_universe(leg, run, stage, lam_val):
    """Generate mda universe object.

    Args:
        leg (str): Bound or free
        run (i): Run number
        stage (str): Restrain, discharge, vanish
        lam_val (str): e.g. "0.000"

    Returns:
        mda.universe: mda universe object
    """
    run_name = dir_paths.get_run_name(run,leg)
    paths = dir_paths.get_dir_paths([run],leg)
    top_file = f"{paths[run_name][stage]['input']}/SYSTEM.top"
    traj_file = f"{paths[run_name][stage]['output']}/lambda-{lam_val}/traj000000001.dcd"
    u = mda.Universe(top_file, traj_file)
    print(f"Opening {traj_file}")

    return u


def read_rest_dict(cfg_file, rest_type):
    """Read config file to extract restraints dict.

    Args:
        cfg_file (str): Path to config file
        rest_type (str): Type of dictionary to read, of multiple_dist, 
        Boresch, or Cart

    Returns:
        rest_dict: Restraints dictionary, as supplied in the config file
    """
    with open(cfg_file,"r") as istream:
        lines = istream.readlines()

    rest_dict = {}
    for l in lines:
        if rest_type == "multiple_dist":
            if l.startswith("distance restraints dictionary"):
                dict_as_list = l.split("=")[1][1:-1] # remove leading space and \n
                rest_dict = literal_eval("".join(dict_as_list))
                break

        if rest_type == "Boresch":
            if l.startswith("boresch restraints dictionary"):
                dict_as_list = l.split("=")[1][1:] # remove leading space 
                if dict_as_list[-1] == "\n":
                    dict_as_list = dict_as_list[:-1] # Remove \n
                rest_dict = literal_eval("".join(dict_as_list))
                break
        
        if rest_type == "Cart":
            if l.startswith("cartesian restraints dictionary"):
                dict_as_list = l.split("=")[1][1:-1] # remove leading space and \n
                rest_dict = literal_eval("".join(dict_as_list))
                break

    return rest_dict


# Functions to get restrained DOF

def get_distance(idx1, idx2, u):
    """Distance between two atoms in Angstrom

    Args:
        idx1 (int): Index of first atom
        idx2 (int): Index of second atom
        u (mda universe): System

    Returns:
        float: Distance in Angstrom
    """
    distance = dist(mda.AtomGroup([u.atoms[idx1]]), mda.AtomGroup([u.atoms[idx2]]), box=u.dimensions)[2][0]
    return distance

def get_angle(idx1, idx2, idx3, u):
    """Angle between three particles in rad.

    Args:
        idx1 (int): Index of first atom
        idx2 (int): Index of second atom
        idx3 (int): Index of third atom
        u (mda universe): System

    Returns:
        float: Angle in rad
    """
    C = u.atoms[idx1].position 
    B = u.atoms[idx2].position 
    A = u.atoms[idx3].position 
    BA = A - B
    BC = C - B
    angle = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
    return angle

def get_dihedral(idx1, idx2, idx3, idx4, u):
    """Get dihedral based on four atom positions

    Args:
        idx1 (int): Index of first atom
        idx2 (int): Index of second atom
        idx3 (int): Index of third atom
        idx4 (int): Index of fourth atom
        u (mda universe): System

    Returns:
        float: Dihedral angle in rad
    """
    positions =[u.atoms[idx].position for idx in [idx1,idx2,idx3,idx4]]
    dihedral = calc_dihedrals(positions[0], positions[1], positions[2], positions[3], box = u.dimensions)
    return dihedral

def get_boresch_dof(anchor_ats, u):
    """Calculate the degrees of freedom defined by the Boresch
    restraints. In addition, get the internal bond angles made 
    by the anchor points in the receptor and ligand, because these
    are important for the stability of the dihedrals. Ordering of
    connection of anchor atoms is r3, r2, r1, l1, l2, l3.

    Args:
        anchor_ats (tuple): Tuple of anchor atom indices of 
        form (r1,r2,r3,l1,l2,l3), where r are anchors on the receptor
        and l on the ligand
        u (mda universe): System

    Returns:
        int, floats: Boresch degrees of freedom
    """
    r1,r2,r3,l1,l2,l3 = anchor_ats
    r = get_distance(r1,l1,u)
    thetaA = get_angle(r2,r1,l1,u)
    thetaB = get_angle(r1,l1,l2,u)
    phiA = get_dihedral(r3,r2,r1,l1,u)
    phiB = get_dihedral(r2,r1,l1,l2,u)
    phiC = get_dihedral(r1,l1,l2,l3,u)
    # Not restrained but distance from coolinearity must be checked
    thetaR = get_angle(r3,r2,r1,u) # Receptor internal angle
    thetaL = get_angle(l1,l2,l3,u) # Ligand internal angle
    return r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL


def get_multiple_dist_dof(multiple_dist_dict, u):
    """Calculate the distances restrained by the multiple distance
    restraints.

    Args:
        multiple_dist_dict (dict): The dictionary specifying the multiple
        distance restraints
        u (mda universe): System

    Returns:
        dict: dictionary of restrained pairs and the distances between them
    """
    multiple_dist_distances = {}
    for pair in multiple_dist_dict:
        multiple_dist_distances[pair] = get_distance(pair[0], pair[1], u)
    return multiple_dist_distances


def track_boresch_dof(anchor_ats, u, percent_traj):
    """Get values, mean, and standard deviation of Boresch
    degrees of freedom and internal angles defined by supplied
    anchor atoms. Also calculate total variance accross all DOF
    , neglecting the internal angles.

    Args:
        anchor_ats (tuple): Anchor atom indices, of form (r1,r2,r3,l1,l2,l3)
        u (mda universe): The system
        percent_traj (float): Percentage of run to average over (25 % would
        result in intial 75 % of trajectory being discarded)

    Returns:
        dict: dictionary of form {"tot_var": tot_var, dof1 :{"mean":mean, "sd":sd,
         "values":[...]}, dof2:{...},...}
    """
    
    r1,r2,r3,l1,l2,l3 = anchor_ats
    n_frames = len(u.trajectory)
    first_frame = round(n_frames - ((percent_traj/100)*n_frames)) # Index of first frame to be used for analysis
    
    dof_dict = {}
    dof_list = ["r","thetaA","thetaB","phiA","phiB","phiC","thetaR","thetaL"]
    # Add sub dictionaries for each Boresch degree of freedom
    for dof in dof_list:
        dof_dict[dof]={}
        dof_dict[dof]["values"]=[]

    for i, frame in enumerate(u.trajectory):
        if i >= first_frame:
            if i == first_frame:
                print(f"First frame no: {i+1}")
            r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL = get_boresch_dof(anchor_ats,u)
            dof_dict["r"]["values"].append(r)
            dof_dict["thetaA"]["values"].append(thetaA)
            dof_dict["thetaB"]["values"].append(thetaB)
            dof_dict["phiA"]["values"].append(phiA)
            dof_dict["phiB"]["values"].append(phiB)
            dof_dict["phiC"]["values"].append(phiC)
            dof_dict["thetaR"]["values"].append(thetaR)
            dof_dict["thetaL"]["values"].append(thetaL)

            if i == n_frames-1:
                print(f"Last frame no: {i+1}")
                dof_dict["tot_var"]=0
                for dof in dof_list:
                    dof_dict[dof]["values"]=np.array(dof_dict[dof]["values"])
                    dof_dict[dof]["mean"]=dof_dict[dof]["values"].mean()
                    # For dihedrals, compute variance based on list of values 
                    # corrected for periodic boundary at pi radians, because there
                    #  is no problem with dihedrals in this region
                    if dof[:3] == "phi":
                        mean = dof_dict[dof]["mean"]

                        # correct variance - fully rigorous
                        corrected_values_sd = []
                        for val in dof_dict[dof]["values"]:
                            dtheta = abs(val - mean)
                            corrected_values_sd.append(min(dtheta, 2*np.pi-dtheta))
                        corrected_values_sd = np.array(corrected_values_sd) 
                        dof_dict[dof]["sd"]=corrected_values_sd.std()

                        # Correct mean (not exact and will fail if very well split above 
                        # and below 2pi)get middle of interval based on current mean
                        print("WARNING: Mean for dihedrals may be incorrect if true mean" \
                              " is near periodic boundary")
                        corrected_values_mean=[]
                        periodic_bound = mean - np.pi
                        if periodic_bound < -np.pi:
                            periodic_bound+=2*np.pi
                        # shift vals from below periodic bound to above
                        for val in dof_dict[dof]["values"]:
                            if val < periodic_bound:
                                corrected_values_mean.append(val+2*np.pi)
                            else:
                                corrected_values_mean.append(val)
                        corrected_values_mean = np.array(corrected_values_mean)
                        mean_corrected = corrected_values_mean.mean()
                        #shift mean back to normal range
                        if mean_corrected > np.pi:
                            dof_dict[dof]["mean"]=mean_corrected-2*np.pi
                        else:
                            dof_dict[dof]["mean"]=mean_corrected
                            
                    else:
                        dof_dict[dof]["sd"]=dof_dict[dof]["values"].std()
                    # Exclude variance of internal angles as these are not restrained
                    if (dof != "thetaR" and dof != "thetaL"):
                        dof_dict["tot_var"]+=dof_dict[dof]["sd"]**2
                    # Assume Gaussian distributions and calculate "equivalent"
                    # force constants for harmonic potentials
                    # so as to reproduce these distributions
                    dof_dict[dof]["k_equiv"]=0.593/(dof_dict[dof]["sd"]**2) # RT at 298 K is 0.593 kcal mol-1
    
    return dof_dict


def track_multiple_dist_dof(multiple_dist_dict, u, percent_traj):
    """Get values, mean, and standard deviation of distances restrained
    using multiple distance restraints.Also calculate total variance accross
    all distances.

    Args:
        multiple_dist_dict (dict): The dictionary specifying the multiple
        distance restraints
        u (mda universe): The system
        percent_traj (float): Percentage of run to average over (25 % would
        result in intial 75 % of trajectory being discarded)

    Returns:
        dict: dictionary of form {"tot_var": tot_var, (ligand_anchor1, protein_anchor1):
        {"mean":mean, "sd":sd, "values":[...]}, (ligand_anchor2, protein_anchor2):{...},...}
    """
    
    n_frames = len(u.trajectory)
    first_frame = round(n_frames - ((percent_traj/100)*n_frames)) # Index of first frame to be used for analysis

    dof_dict = {}
    # Add sub dictionaries for each distance
    for pair in multiple_dist_dict:
        dof_dict[pair]={}
        dof_dict[pair]["values"]=[]

    for i, frame in enumerate(u.trajectory):
        if i >= first_frame:
            if i == first_frame:
                print(f"First frame no: {i+1}") # Convert index to number
            multiple_dist_dists = get_multiple_dist_dof(multiple_dist_dict, u)
            for pair in multiple_dist_dists:
                dof_dict[pair]["values"].append(multiple_dist_dists[pair])

            if i == n_frames-1:
                print(f"Last frame no: {i+1}")
                dof_dict["tot_var"]=0
                for pair in dof_dict:
                    if pair != "tot_var": # Check pair is actually pair
                        dof_dict[pair]["values"]=np.array(dof_dict[pair]["values"])
                        dof_dict[pair]["mean"]=dof_dict[pair]["values"].mean()
                        dof_dict[pair]["sd"]=dof_dict[pair]["values"].std()
                        dof_dict["tot_var"]+=dof_dict[pair]["sd"]**2
                        # Assume Gaussian distributions and calculate "equivalent"
                        # force constants for harmonic potentials
                        # so as to reproduce these distributions
                        dof_dict[pair]["k_equiv"]=0.593/(dof_dict[pair]["sd"]**2) # RT at 298 K is 0.593 kcal mol-1
    
    return dof_dict


# Functions required for tracking Cartesian degrees of freeedom

def get_coord_sys(idx1, idx2, idx3, u):
    # Origin
    orig = u.atoms[idx1].position
    # Get x axis
    X_unnorm = u.atoms[idx2].position - u.atoms[idx1].position
    X = X_unnorm/norm(X_unnorm)
    # Get y axis plane - not orthogonal to X
    Y_plane = u.atoms[idx3].position - u.atoms[idx2].position # Lies in x-y plane but not orthogonal to x
    # Z axis is perpendicular to plane made by X and Y_plane
    Z_unnorm = np.cross(X, Y_plane)
    Z = Z_unnorm/norm(Z_unnorm)
    # Compute Y to be orthogonal to X and Z
    Y_unnorm = -np.cross(X, Z) # Negative because need to use left-hand rule. Should be normalised already, but normalise anayway.
    Y = Y_unnorm/norm(Y_unnorm)

    return orig, X, Y, Z   

class CoordSys:
    def __init__(self, idx1, idx2, idx3, u):
        self.idx1, self.idx2, self.idx3 = idx1, idx2, idx3
        self.orig, self.X, self.Y, self.Z = get_coord_sys(idx1, idx2, idx3, u)
        self.offset_coords = None
        self.already_offset = False

    def get_point_pos(self, point):
        """Get the position of a point in the coordinate system"""
        # Get vector to point based on coordinate system origin
        P = point - self.orig
        x = np.dot(self.X, P)
        y = np.dot(self.Y, P)
        z = np.dot(self.Z, P)
        return(x, y, z)

    def get_coord_sys_pos(self, ref_coord_sys):
        """Get position of reference coordinate system in frame of
        reference of current coordinate system

        Args:
            ref_coord_sys (CoordSys): The reference coordinate system whose
            position is to be described by the current coordinate system
        """
        positions = {"X":{}, "Y":{}, "Z":{}}
        positions["X"]["local_x"] = np.dot(ref_coord_sys.X, self.X) 
        positions["X"]["local_y"] = np.dot(ref_coord_sys.X, self.Y) 
        positions["X"]["local_z"] = np.dot(ref_coord_sys.X, self.Z) 
        positions["Y"]["local_x"] = np.dot(ref_coord_sys.Y, self.X)
        positions["Y"]["local_y"] = np.dot(ref_coord_sys.Y, self.Y)
        positions["Y"]["local_z"] = np.dot(ref_coord_sys.Y, self.Z)
        positions["Z"]["local_x"] = np.dot(ref_coord_sys.Z, self.X)
        positions["Z"]["local_y"] = np.dot(ref_coord_sys.Z, self.Y)
        positions["Z"]["local_z"] = np.dot(ref_coord_sys.Z, self.Z)
        return positions
    
    def offset(self, ref_frame_rot):
        """Offset the coordinate system to superimpose on a reference
        coordinate system.

        Args:
            ref_frame_rot (dict): Dictionary giving reference frame rotation in same form as in config file.
            u (mda.universe): Universe
        """
        if not self.already_offset:
            temp_X = self.X*ref_frame_rot["xl_ref"]["xl_ref_xl"] + self.Y*ref_frame_rot["xl_ref"]["xl_ref_yl"] + self.Z*ref_frame_rot["xl_ref"]["xl_ref_zl"]
            temp_Y = self.X*ref_frame_rot["yl_ref"]["yl_ref_xl"] + self.Y*ref_frame_rot["yl_ref"]["yl_ref_yl"] + self.Z*ref_frame_rot["yl_ref"]["yl_ref_zl"]
            temp_Z = self.X*ref_frame_rot["zl_ref"]["zl_ref_xl"] + self.Y*ref_frame_rot["zl_ref"]["zl_ref_yl"] + self.Z*ref_frame_rot["zl_ref"]["zl_ref_zl"]
            self.X, self.Y, self.Z = temp_X, temp_Y, temp_Z

            #print(self.X, self.Y, self.Z)
            #print(self.X, self.Y, self.Z)
            #try:
                #self.check_valid(self.X, self.Y, self.Z)
            #except Exception as e:
                #print(e.with_traceback)
                #print(offset_coords)
            self.check_valid(self.X, self.Y, self.Z)
            self.ref_frame_rot = ref_frame_rot
            self.already_offset = True
        else:
            raise Exception("Axis has already been offset") 


    def update(self, u):
        self.orig, self.X, self.Y, self.Z = get_coord_sys(self.idx1, self.idx2, self.idx3, u)
        if self.already_offset:
            temp_X = self.X*self.ref_frame_rot["xl_ref"]["xl_ref_xl"] + self.Y*self.ref_frame_rot["xl_ref"]["xl_ref_yl"] + self.Z*self.ref_frame_rot["xl_ref"]["xl_ref_zl"]
            temp_Y = self.X*self.ref_frame_rot["yl_ref"]["yl_ref_xl"] + self.Y*self.ref_frame_rot["yl_ref"]["yl_ref_yl"] + self.Z*self.ref_frame_rot["yl_ref"]["yl_ref_zl"]
            temp_Z = self.X*self.ref_frame_rot["zl_ref"]["zl_ref_xl"] + self.Y*self.ref_frame_rot["zl_ref"]["zl_ref_yl"] + self.Z*self.ref_frame_rot["zl_ref"]["zl_ref_zl"]
            self.X, self.Y, self.Z = temp_X, temp_Y, temp_Z
        self.check_valid(self.X, self.Y, self.Z)


    def get_euler_angles(self, ref_coord_sys, u):
        """Get the Euler angles (factored as RzRyRx, see implementation
        in Sire) describing the rotation of the reference coordinate
        system with respect to the current system.

        Args:
            ref_coord_sys (CoordSys): Reference coordinate system
            u (mda.universe): Universe

        Returns:
            floats: Euler angles phi, theta, psi
        """
        positions = self.get_coord_sys_pos(ref_coord_sys)
        phi = np.arctan2(positions["Y"]["local_x"], positions["X"]["local_x"]) # [-pi, pi]
        theta = np.arcsin(- positions["Z"]["local_x"]) # [0, pi]
        psi = np.arctan2(positions["Z"]["local_y"], positions["Z"]["local_z"]) # [-pi, pi]

        return phi, theta, psi

    @staticmethod
    def check_valid(X, Y, Z):
        # Norms == 1
        assert(round(norm(X), 4) == 1.0)
        assert(round(norm(Y), 4) == 1.0)
        assert(round(norm(Y), 4) == 1.0)
        # Orthogonal
        assert(round(np.dot(X, Y), 4) == 0)
        assert(round(np.dot(X, Z), 4) == 0)
        assert(round(np.dot(Y, Z), 4) == 0)

        return True
        

def get_cart_dof(lig_sys, recept_sys, u):
    xr_l1, yr_l1, zr_l1 = recept_sys.get_point_pos(lig_sys.orig)
    phi, theta, psi = recept_sys.get_euler_angles(lig_sys, u)
    return xr_l1, yr_l1, zr_l1, phi, theta, psi


def track_cartesian_dof(anchor_ats, ref_frame_rot, u, percent_traj):
    """Get values, mean, and standard deviation of Cartesian
    degrees of freedom and internal angles defined by supplied
    anchor atoms. Also calculate total variance accross all DOF
    , neglecting the internal angles.

    Args:
        anchor_ats (tuple): Anchor atom indices, of form (r1,r2,r3,l1,l2,l3)
        ref_frame_rot (dict): Reference frame rotation dictionary, in same form as in config file
        u (mda universe): The system
        percent_traj (float): Percentage of run to average over (25 % would
        result in intial 75 % of trajectory being discarded)

    Returns:
        dict: dictionary of form {"tot_var": tot_var, dof1 :{"mean":mean, "sd":sd,
         "values":[...]}, dof2:{...},...}
    """
    
    l1, l2, l3, r1, r2, r3 = anchor_ats
    n_frames = len(u.trajectory)
    first_frame = round(n_frames - ((percent_traj/100)*n_frames)) # Index of first frame to be used for analysis
    
    dof_dict = {}
    dof_list = ["xr_l1","yr_l1","zr_l1","phi","theta","psi","thetaL","thetaR"]
    # Add sub dictionaries for each Boresch degree of freedom
    for dof in dof_list:
        dof_dict[dof]={}
        dof_dict[dof]["values"]=[]

    # Get Coordinate systems and set offset for ligand coordinate system
    recept_sys = CoordSys(r1, r2, r3, u)
    lig_sys = CoordSys(l1, l2, l3, u)
    lig_sys.offset(ref_frame_rot)

    # Populate these dictionaries with values from trajectory
    n_frames = len(u.trajectory)

    for i, frame in enumerate(u.trajectory):
        if i >= first_frame:
            if i == first_frame:
                print(f"First frame no: {i+1}")
            lig_sys.update(u)
            recept_sys.update(u)
            xr_l1, yr_l1, zr_l1, phi, theta, psi = get_cart_dof(lig_sys, recept_sys, u)
            thetaL = get_angle(r1, r2, r3, u)
            thetaR = get_angle(l1, l2, l3, u)
            dof_dict["xr_l1"]["values"].append(xr_l1)
            dof_dict["yr_l1"]["values"].append(yr_l1)
            dof_dict["zr_l1"]["values"].append(zr_l1)
            dof_dict["phi"]["values"].append(phi)
            dof_dict["theta"]["values"].append(theta)
            dof_dict["psi"]["values"].append(psi)
            dof_dict["thetaL"]["values"].append(thetaL)
            dof_dict["thetaR"]["values"].append(thetaR)

            if i == n_frames-1:
                print(f"Last frame no: {i+1}")
                dof_dict["tot_var"]=0
                for dof in dof_list:
                    dof_dict[dof]["values"]=np.array(dof_dict[dof]["values"])
                    dof_dict[dof]["mean"]=dof_dict[dof]["values"].mean()
                    dof_dict[dof]["sd"]=dof_dict[dof]["values"].std()
                    # Exclude variance of internal angles as these are not restrained
                    if (dof != "thetaR" and dof != "thetaL"):
                        dof_dict["tot_var"]+=dof_dict[dof]["sd"]**2
                    # Assume Gaussian distributions and calculate "equivalent"
                    # force constants for harmonic potentials
                    # so as to reproduce these distributions
                    dof_dict[dof]["k_equiv"]=0.593/(dof_dict[dof]["sd"]**2) # RT at 298 K is 0.593 kcal mol-1
    
    return dof_dict


def get_dof_dicts(leg, runs, stage, lam_val, percent_traj, dof_type):
    """Get dof_dicts for given stage and lambda value
    for all supplied runs

    Args:
        leg (str): bound
        runs (list): Run numbers (ints)
        stage ([type]): restrain, discharge, or vanish
        lam_val (str): Window of interest
        percent_traj (float): Percentage of run to average over (25 % would
        result in intial 75 % of trajectory being discarded)
        dof_type (str): Boresch, multiple_dist, Cart

    Returns:
        dict: dict of dof_dicts, with run names as keys
    """
    if type(lam_val) == float:
        lam_val = f"{lam_val:.3f}" # Account for input of float instead of string
    dof_dicts = {}
    paths = dir_paths.get_dir_paths(runs, leg)

    for run in runs:
        run_name = dir_paths.get_run_name(run,leg)
        dof_dicts[run_name] = {}
        u = get_mda_universe(leg, run, stage, lam_val)
        cfg_path = f'{paths[run_name][stage]["input"]}/sim.cfg'
        dof_dict = {}
        if dof_type == "Boresch":
            rest_dict = read_rest_dict(cfg_path, rest_type=dof_type)
            anchor_ats = tuple([x for x in rest_dict["anchor_points"].values()])
            dof_dict = track_boresch_dof(anchor_ats, u, percent_traj)
        if dof_type == "multiple_dist":
            rest_dict = read_rest_dict(cfg_path, rest_type=dof_type)
            dof_dict = track_multiple_dist_dof(rest_dict, u, percent_traj)
        elif dof_type == "Cart":
            rest_dict = read_rest_dict(cfg_path, rest_type=dof_type)
            anchor_ats = tuple([x for x in rest_dict["anchor_points"].values()])
            ref_frame_rot = rest_dict["reference_frame_rotation"]
            dof_dict = track_cartesian_dof(anchor_ats, ref_frame_rot, u, percent_traj)
        dof_dicts[run_name]=dof_dict

    return dof_dicts
    

def plot_dof_hists(leg, runs, stage, lam_val, percent_traj, selected_dof_list, dof_type):
    """Plot histograms of selected degrees of freedom over specified
    final percentage of trajectory for supplied runs and lambda window.

    Args:
        leg (str): bound
        runs (list): Run numbers (ints)
        stage ([type]): restrain, discharge, or vanish
        lam_val (str): Window of interest
        percent_traj (float): Percentage of run to average over (25 % would
        result in intial 75 % of trajectory being discarded)
        selected_dof_list (list): For Boresch restraints, subset of ["r","thetaA",
        "thetaB","phiA","phiB","phiC","thetaR","thetaL"] to be plotted. For multiple
        distance restraints, this defaults to all pairs if the supplied list is empty.
        dof_type (str): Boresch, multiple_dist, Cartesian
    """
    print("###############################################################################################")
    print(f"Plotting histograms for {dof_type} DOF for {leg} {stage} lambda = {lam_val} and final {percent_traj} % of traj")

    dof_dicts = get_dof_dicts(leg, runs, stage, lam_val, percent_traj, dof_type)
    if dof_type == "multiple_dist" and selected_dof_list == []: # Select all pairs
        run_names = [run_name for run_name in dof_dicts]
        selected_dof_list = [pair for pair in dof_dicts[run_names[0]] if pair != "tot_var"]
    no_dof = len(selected_dof_list)

    #fig, axs = plt.subplots(ceil(no_dof/6), 6, figsize=(4*6,4*ceil(no_dof/6)), dpi=500)
    fig, axs = plt.subplots(ceil(no_dof/3), 3, figsize=(4*3,2*ceil(no_dof/3)), dpi=800)
    #colours =  ['#00429d', '#3a8bbb', '#ffb59a', '#ff6b95', '#93003a'] # Will cause issues for more than 5 runs
    colours =  ['#000000', '#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00'] # Will cause issues for more than 10 runs
    axs = axs.flatten()

    # Dictionary to translate the Boresch DOF to latex for figure labels
    boresch_to_latex = {"r":r"$r_{Aa}$","thetaA":r"$\theta_A$","thetaB":r"$\theta_B$","phiA":"$\phi_A$","phiB":"$\phi_B$",
        "phiC":"$\phi_C$","thetaR":r"$\theta_R$","thetaL":r"$\theta_L$"}

    for j, run in enumerate(runs):
        run_name = dir_paths.get_run_name(run,leg)
        for i, dof in enumerate(selected_dof_list):
            ax = axs[i]
            values = dof_dicts[run_name][dof]["values"]
            mean = dof_dicts[run_name][dof]["mean"]
            sd = dof_dicts[run_name][dof]["sd"]
            ax.hist(values,label = f"{run_name}", color=colours[j], edgecolor='k')
            ax.axvline(mean, linestyle = "dashed", color=colours[j], linewidth=2, label=f"Mean: {mean:.2f}\nSD: {sd:.2f}")
            if dof == "r" or dof[1] == "r":
                ax.set_ylabel(f"{boresch_to_latex[dof]} " + r"($\mathrm{\AA}$)")
            elif type(dof) == tuple:
                ax.set_ylabel(f"Dist between indices {dof[0]} and {dof[1]}" + r" $\textrm{\AA}$")
            else:
                ax.set_ylabel(f"{boresch_to_latex[dof]} (rad)")
            ax.set_ylabel("Counts")
            ax.legend(loc=(1.04,0))

    fig.tight_layout()
    mkdir_if_required("analysis")
    mkdir_if_required(f"analysis/{dof_type}_dof")
    fig.savefig(f"analysis/{dof_type}_dof/{leg}_{stage}_{lam_val:.3f}_{dof_type}_dof_hists.png")


def plot_dof_vals(leg, runs, stage, lam_val, percent_traj, selected_dof_list, dof_type):
    """Plot values of selected degrees of freedom over specified
    final percentage of trajectory for supplied runs and lambda window.

    Args:
        leg (str): bound
        runs (list): [
        runs (list): Run numbers (ints)
        stage ([type]): restrain, discharge, or vanish
        lam_val (str): Window of interest
        percent_traj (float): Percentage of run to average over (25 % would
        result in intial 75 % of trajectory being discarded)
        selected_dof_list (list): Subset of ["r","thetaA","thetaB","phiA","phiB",
        "phiC","thetaR","thetaL"]
        dof_type (str): Boresch, multiple_dist, Cart
    """
    print("###############################################################################################")
    print(f"Plotting values of Boresch DOF for {leg} {stage} lambda = {lam_val} and final {percent_traj} % of traj")

    dof_dicts = get_dof_dicts(leg, runs, stage, lam_val, percent_traj, dof_type)
    if dof_type == "multiple_dist" and selected_dof_list == []: # Select all pairs
        run_names = [run_name for run_name in dof_dicts]
        selected_dof_list = [pair for pair in dof_dicts[run_names[0]] if pair != "tot_var"]
    no_dof = len(selected_dof_list)

    fig, axs = plt.subplots(ceil(no_dof/3), 3, figsize=(4*3,2*ceil(no_dof/3)), dpi=800)
    #colours =  ['#00429d', '#3a8bbb', '#ffb59a', '#ff6b95', '#93003a'] # Will cause issues for more than 5 runs
    colours =  ['#000000', '#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00'] # Will cause issues for more than 10 runs
    axs = axs.flatten()

    # Dictionary to translate the Boresch DOF to latex for figure labels
    boresch_to_latex = {"r":r"$r_{Aa}$","thetaA":r"$\theta_A$","thetaB":r"$\theta_B$","phiA":"$\phi_A$","phiB":"$\phi_B$",
        "phiC":"$\phi_C$","thetaR":r"$\theta_R$","thetaL":r"$\theta_L$"}


    for j, run in enumerate(runs):
        run_name = dir_paths.get_run_name(run,leg)
        for i, dof in enumerate(selected_dof_list):
            ax = axs[i]
            values = dof_dicts[run_name][dof]["values"]
            mean = dof_dicts[run_name][dof]["mean"]
            sd = dof_dicts[run_name][dof]["sd"]
            ax.plot([x for x in range(len(values))], values, label = f"{run_name}", color=colours[j])
            ax.axhline(mean, linestyle = "dashed", color=colours[j], linewidth=2, label=f"Mean: {mean:.2f}\nSD: {sd:.2f}")
            if dof == "r" or dof[1] == "r":
                ax.set_ylabel(f"{boresch_to_latex[dof]} " + r"($\mathrm{\AA}$)")
            elif type(dof) == tuple:
                ax.set_ylabel(f"Dist between indices {dof[0]} and {dof[1]}" + r" $\textrm{\AA}$")
            else:
                ax.set_ylabel(f"{boresch_to_latex[dof]} (rad)")
            ax.set_xlabel("Frame No")
            ax.legend(loc=(1.04,0))
#           handles, labels = ax.get_legend_handles_labels()
#           lg = fig.legend(handles, labels, bbox_to_anchor=(1.10, 0.5))
                #ax.legend().set_visible(False)

    fig.tight_layout()
    mkdir_if_required("analysis")
    mkdir_if_required(f"analysis/{dof_type}_dof")
    fig.savefig(f"analysis/{dof_type}_dof/{leg}_{stage}_{lam_val:.3f}_{dof_type}_dof_vals.png")
#                bbox_extra_artists=(lg,), 
#                bbox_inches='tight')


# Just seems to return whitespace (but no errors) if functions modified to return figures
#def plot_dof(leg, runs, stage, lam_val, percent_traj, selected_dof_list):
#    fig = plt.figure(figsize=(4*len(selected_dof_list), 8))
#    subfigs = fig.subfigures(1, 2)
#    subfigs[0] = plot_dof_hists(leg, runs, stage, lam_val, percent_traj, selected_dof_list)
#    subfigs[1] = plot_dof_vals(leg, runs, stage, lam_val, percent_traj, selected_dof_list, False)
#    fig.savefig(f"analysis/{leg}_{stage}_{lam_val}_boresch_dof.png")
