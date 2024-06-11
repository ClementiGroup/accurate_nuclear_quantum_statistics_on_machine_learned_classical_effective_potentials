---
layout: default
---

# Training data generation

 Dataset generation

For a given temperature, the dataset generation involves two steps.
First, PIMD simulation for the system of interest is performed. 
The input files and SLURM ubmission script required to run i-pi simulation can be found in folder `example_setup`.

Second, the raw data are processed:  the i-pi output is combined into a single hdf5 dataset (parameters are  given in `config.yaml`). Then, total forces are calculated (`config_postprocess.yaml`), followed by force optimisaziton and substruction of the prior (`config_portprocess.yaml`)

Since the compute and  memory requirements are different for each of these steps,the postprocessing is done as a sequence of jobs, with the correspondign SLURM submision files `slurm_gen_data.job`, `slurm_postprocess_data.job`, `slurm_optimize.job`

It is assumed, that PIMD simulation are located in folders with names corresponding to temperatures (eg. `300`)


The resulting dataset `h2o_256.h5` can be dounloaded from Zenodo
#TODO: add links to Zenodo repository and instructions on which file to use 