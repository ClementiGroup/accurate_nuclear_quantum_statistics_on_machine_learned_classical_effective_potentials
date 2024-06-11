---
layout: default
---

# Molecular dynamics simulation using the trained model

The folder `example_setup` contains all the files, required for the simulation

`h2o.xyz` : starting configuration for the simulation
`input.xml` : i-pi configuration file.
`run-ase.yaml` : configuration of the model
`submit.job` : SLURM submission script

The setup is given for 100 K. To perform a simulations at a given temperature and/or a different model, the following adjustments need to be made: 

1) `input.xml`:  Set simulation temperature (tag `<temperature>`) and  initial temperature for velosity distribution (tag `<velocities>`) to the target temperature (K)
2) `input.xml` : Set the prior rescaling parameter  in line 27 (` <force forcefield="pes" weight="0.16666666666666666" />`) to $\frac{T}{T_0}$, where $T_0=600 \text{K}$.
3) `run-ase.yaml`:  Set path to the checkpoint file of interest (`PATH2CKPT`)
4) `run-ase.yaml`: Populate field `PATH2PIGS` with the value of `$TEPIGS_PATH` environment variable
5) `submit.job`: Update to reflect the computational resources available
6)  In all the files, replace string `socket_address` with a unique socket address 


