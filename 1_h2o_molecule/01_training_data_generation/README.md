---
layout: default
---

# Training data generation

 Dataset generation

For a given temperature, the dataset generation involves two steps.
First, PIMD simulation for the system of interest is performed. 
An example initialization xyz file, i-pi xml configuration file and slurm submission script required to run i-pi simulation can be found in folder `example_setup`.

Second, the raw data are processed:  the i-pi output is combined into a single hdf5 dataset (parameters are  given in `config.yaml`). Then, total forces are calculated (`config_postprocess.yaml`), followed by force optimisaziton and substruction of the prior (`config_portprocess.yaml`)

The entire prostprocessing step is handled by `slurm.job`
It is assumed, that PIMD simulation are located in folders with names corresponding to temperatures in the following format {T:>04d} (eg. `0100`, `0300`, `0600`)


The resulting dataset `h2o.h5` can be dounloaded from Zenodo
#TODO: add links to Zenodo repository and instructions on which file to use 
