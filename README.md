# Alchemical Simulations with BioSimSpace

## Aimed at

Anyone interested in learning to perform alchemical free energy calculations using the CCP-BioSim
[BioSimSpace](https://github.com/michellab/BioSimSpace) Python environment for easy setup, running, management and analysis of biomolecular simulations.

## Requirements

Knowledge of Python, e.g. as gained from the
<a href="https://training.ccpbiosim.ac.uk/workshops/hub/spawn?profile=python-workshop" target="_blank">Python for Biomolecular Modellers</a> workshop.

Basic knowledge of atomistic simulations.

## Abstract

Alchemical free energy calculations can be used to efficiently compute binding free energies between a ligand and a protein or hydration free energies of a small molecule. In the last few years, the use of such methods has gained momentum not only within academia but also within the pharmaceutical industry. In order to run alchemical free energy simulations, a series of molecular dynamics simulations need to be carried out. During this workshop you will learn how to set up, run, and analyse both relative and absolute binding free energy calculations with BioSimSpace.

## Training Material: 

The workshop consists of a series of Jupyter notebooks. These are available on the
<a href="https://workshop.biosimspace.org" target="_blank">workshop jupyter server</a>
and can be downloaded from the <a href="https://github.com/michellab/bssccpbiosim2022" target="_blank">GitHub repository</a>.

Once you have started the server, navigate to the `workshops/` directory and you will find the
notebooks there. These training materials will teach you more about BioSimSpace, including how to write your own BioSimSpace code. The material is split into three parts.

### Part 1: Introduction to Alchemistry with BioSimSpace

* [alchemical_setup.ipynb](introduction_to_alchemistry/alchemical_introduction.ipynb) : This notebook shows you how to use BioSimSpace to set up alchemical free energy simulations that can be run with [SOMD](https://siremol.org/tutorials/somd) or [GROMACS](http://www.gromacs.org).

### Part 2: Relative Binding Free Energies with BioSimSpace

* [setup_rbfe.ipynb](relative_binding_free_energies/setup_rbfe.ipynb) : This notebook shows how to prepare and deploy a network of relative binding free energy calculations. 

* [analyse_rbfe.ipynb](relative_binding_free_energies/analysis_rbfe.ipynb) : This notebook shows how to analyse a network of relative binding free energies to generate estimates of protein-ligand binding affinities. 

### Part 3: Absolute Binding Free Energies with BioSimSpace

* [setup_abfe.ipynb](absolute_binding_free_energies/setup_abfe.ipynb) : This notebook shows how to prepare an absolute binding free energy calculation for a protein-ligand complex using BioSimSpace.

* [analyse_rbfe.ipynb](absolute_binding_free_energies/analyse_abfe.ipynb) :This notebook shows how to analyse an absolute binding free energy calculation of a protein-ligand complex using BioSimSpace.

