---
layout: default
---
# Model training 

To perform the training, download the dataset from [Zenodo](https://doi.org/10.5281/zenodo.12684727). Unpack the `CG_quantum_statistics.zip` archive and copy files from `CG_quantum_statistics/0_morse_potential` to  `accurate_nuclear_quantum_statistics_on_machine_learned_classical_effective_potentials/0_morse_potential/01_training_data_generation/`

The folder `example_setup` contains files required for training:  

`config.yaml` : parameters of the model. Path to the position and force file should be adjusted accordingly.  

`slurm.job` : SLURM submission file


