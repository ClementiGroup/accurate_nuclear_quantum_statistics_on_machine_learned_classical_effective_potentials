---
layout: default
---

# Dataset generation

For a given temperature, the dataset generation involves two steps.

First, PIMD simulation for a single particle in Morse potential is performed. Required initialization xyz file, i-pi xml configuration file and slurm submission script required to run i-pi simulation can be found in folder `example_setup`.

Second,  the i-pi output is combined into a dataset that is then used in training. The processing is handled by `slurm.job`, the processing parameters are
given in `config_template.yaml`  
It is assumed, that PIMD simulation are located in folders with names corresponding to temperatures (eg. `100`, `300`, `600`)


The resulting dataset can be dounloaded from Zenodo
#TODO: add links to Zenodo repository and instructions on which file to use 
