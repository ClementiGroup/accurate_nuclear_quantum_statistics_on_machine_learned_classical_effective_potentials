---
layout: default
---

# Classical MD simulation with I-PI

`example_setup` folder contains files, required for classical MD simulation of a Zundel cation with i-pi:

 - `h5o2.dms4B.coeff.com.dat` and `h5o2.pes4B.coeff.dat` :  parameters of the HBB dipole and potential energy surface of Zundel cation as described in  [J. Chem. Phys. 122, 044308 (2005)](https://doi.org/10.1063/1.1834500).

 - `init.xyz` : Initial system configuration
 - `input.xml`: I-PI configuration file
 - `run-driver.sh`: call to a system-specific driver
 - `submit.sh` : SLURM submission file

 

## Adjustments for different temperature 

The provided example can be used to run the simulation at 100 K. To run simulation at a different temperature, the following adjustments should be made: 

Parameters to change in input files for running the simulation:
`input.xml`:  line  21   : Set correct temperature instead of 100 
            line  17   : Use temperature as a suffix for socket address (replace 100 with correct temperature)

`submit.sh`: lines 15-16      : Set correct suffix for the file to be removed  (replace 100 with correct temperature)

`run-driver.sh`:
           line 1     : Set the same socket address as in line 17 of `input.xml`
