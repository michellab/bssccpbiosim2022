# Function to save data and create required directory 

import os
import matplotlib.pyplot as plt

def mkdir_if_required(name):
    """Create directory if it does not exist

    Args:
        name (str): Name of directory 
    """
    if not os.path.exists(name):
        os.makedirs(name)
