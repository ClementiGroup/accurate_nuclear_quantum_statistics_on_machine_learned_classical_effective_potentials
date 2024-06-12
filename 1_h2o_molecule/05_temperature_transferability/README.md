---
---
# Temperature transferability in the ground state

A model trained on data, generated at some temperature $T_0$ can be used to run simulations at any temperature $T \< T_0$, providing that the system of interest is 
in the ground vibrational state at $T_0$.


The folder `example_setup` contains all the files, required to reproduce the results
presented in the paper 

`h2o.xyz` : starting configuration for the simulation
`input.xml` : i-pi configuration file.
`run-ase.yaml` : configuration of the model
`submit.job` : SLURM submission script

The setup is given for the case, where the model was trained with data produced at T'= 600 K,  and the simulation is done at T=100 K. 


To perform a simulations at a given temperature (T), using a model trained with data generated at T', the  the following adjustments need to be made: 

1) `input.xml`:  Set simulation temperature (tag `<temperature>`) and  initial temperature for velosity distribution (tag `<velocities>`) to the target temperature  T(K)
2) `input.xml` : Set the prior rescaling parameter  in line 27 (` <force forcefield="pes" weight="0.16666666666666666" />`) to $\frac{T}{T_0}$, where $T_0=600 \text{K}$.
3) `input.xml` : Set the prior rescaling parameter  in line 26 (` <force forcefield="fes" weight="0.16666666666666666" />`) to $\frac{T}{T'}$
4) `run-ase.yaml`:  Set path to the checkpoint file of interest (`PATH2CKPT`)
5) `run-ase.yaml`: Populate field `PATH2PIGS` with the value of `$PATH2PIGS` environment variable
6) `submit.job`: Update to reflect the computational resources available
7)  In all the files, replace string `socket_address` with a unique socket address 


