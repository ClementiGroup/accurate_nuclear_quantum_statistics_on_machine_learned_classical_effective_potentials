---
layout: default
---

# Training data generation

 Dataset generation

For a given temperature, the dataset generation involves two steps.
First, PIMD simulation for the system of interest is performed. 
An example initialization xyz file, i-pi xml configuration file and slurm submission script required to run i-pi simulation can be found in folder `example_setup`. 
Files `h5o2.dms4B.coeff.com.dat` and `h5o2.pes4B.coeff.dat` contain parameters of the HBB dipole and potential energy surface of Zundel cation as described in  [J. Chem. Phys. 122, 044308 (2005)](https://doi.org/10.1063/1.1834500).


Second, the raw data are processed:  the i-pi output is combined into a single hdf5 dataset (parameters are  given in `config.yaml`). Then, total forces are calculated (`config_postprocess.yaml`), followed by force optimisaziton and substruction of the prior (`config_portprocess.yaml`)

The entire prostprocessing step is handled by `slurm.job`
It is assumed, that PIMD simulation are located in folders with names corresponding to temperatures (eg. `100`, `300`)


The resulting dataset `h5o2+.h5` can be dounloaded from Zenodo
#TODO: add links to Zenodo repository and instructions on which file to use 
