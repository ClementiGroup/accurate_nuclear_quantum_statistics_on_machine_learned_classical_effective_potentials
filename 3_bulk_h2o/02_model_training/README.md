---
layout: default
---
# Model training 
For training, the environmental variable `$PATH2PIGS` should be set to the location of the `te-pigs` 

To perform the training, download the dataset from [Zenodo](https://doi.org/10.5281/zenodo.12684727). Unpack the `CG_quantum_statistics.zip` archive and copy `CG_quantum_statistics/3_bulk_h2o/h2o_256.h5` to  `accurate_nuclear_quantum_statistics_on_machine_learned_classical_effective_potentials/3_bulk_h2o/01_training_data_generation/`
The folder `example_setup` contains files required for training:  

`mace.yaml` : parameters of the model. These parameters remain  the same for all the models reported.

`partition_settings_global.yml` : parameters of test train split. Temperature and seeds should be adjusted as necessary

`slurm.job` : SLURM submission file

`trainer.yaml` : parameters of the trainings. Temperature, seeds and forces need to be adjusted accordingly. The forces
should be computed with prior scaling that matches the temperature: 
 $\alpha = \frac{T}{T_0}$ , where $T_0=600 \text{K}$. Thus, for 300 K the dataset  'total_forces_projected_scaled_classical_prior_scaling_5.000e-01_l2reg_5.000e-02'  should be used.

For initialization, the following seeds were used: 
1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239
