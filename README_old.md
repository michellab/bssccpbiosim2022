# Advanced Simulations with BioSimSpace

## Aimed at

Anyone interested in learning how to use the new CCP-BioSim
[BioSimSpace](https://github.com/michellab/BioSimSpace) Python environment for easy setup, running, management and analysis of advanced biomolecular simulations.

## Requirements

Knowledge of Python, e.g. as gained from the
<a href="https://ccpbiosim.github.io/python_and_data" target="_blank">Python for Biomolecular Modellers</a> workshop.

Knowledge of BioSimSpace, e.g. as gained from the
<a href="https://ccpbiosim.github.io/biosimspace_workshop" target="_blank">Introduction to BioSimSpace</a> workshop.

## Abstract

Alchemical free energy calculations can be used to efficiently compute binding free energies between a ligand and a protein or hydration free energies of a small molecule. In the last few years, the use of such methods has gained momentum not only within academia but also within the pharmaceutical industry. In order to run alchemical free energy simulations, a series of molecular dynamics simulations need to be carried out. In the first part of this workshop you will learn how to set up, run, and analyse, binding free energy calculations with BioSimSpace.

Metadynamics is an advanced simulation method that can be used to bypass large free energy barriers allowing the estimation of free energies for processes that would normally be impossible to observe. The method relies on the concept of "collective variables", which are used to describe the pathways between basins on the free energy surface, e.g. the distance between the centre of mass of a ligand and the binding site in a protein. In the second part of this workshop you will lean how to use BioSimSpace to set up metadynamics simulations for simple collective variables that can be run using [GROMACS](http://www.gromacs.org) and [PLUMED](https://www.plumed.org).

## Training Material

The workshop consists of a series of Jupyter notebooks. These are available on the
<a href="https://notebook.biosimspace.org" target="_blank">workshop jupyter server</a>
and can be downloaded from the <a href="https://github.com/ccpbiosim/biosimspace-advanced-simulation" target="_blank">GitHub repository</a>.

Once you have started the server, navigate to the `workshops/advanced` directory and you will find the
notebooks there. These training materials will teach you more about BioSimSpace, including how to write your own BioSimSpace code. The material is split into two parts.

### Part 1: Alchemical Free Energy Setup

* [alchemical_setup.ipynb](alchemistry/html/Alchemical_setup.html) : This notebook shows you how to use BioSimSpace to set up alchemical free energy simulations that can be run with [SOMD](https://siremol.org/tutorials/somd) or [GROMACS](http://www.gromacs.org).

### Part 2: Metadynamics

* [metadynamics.ipynb](metadynamics/html/metadynamics.html) : This notebook shows you how to use BioSimSpace to setup, run, and analyse metadynamics simulations using [GROMACS](http://www.gromacs.org) and [PLUMED](https://www.plumed.org).
