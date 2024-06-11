---
layout: default
---

# Classical MD simulation with I-PI


`example_setup` folder contains files, required for classical MD simulation of a single water molecule with i-pi:

 - `init.xyz` : Initial system configuration
 - `input.xml`: I-PI configuration file
 - `run-driver.sh`: call to a system-specific driver
 - `submit.sh` : SLURM submission file

 

Parameters to change in input files for running the simultation:
`input.xml`:  line  21   : Set correct temperature instead of 100 
            line  17   : Use temperature as a suffix for socket address (replace 100 with correct temperature)

`submit.sh`: lines 14-15      : Set correct suffix for the file to be removed  (replace 100 with correct temperature)

`run-driver.sh`:
           line 1     : Set the same socket address as in line 17 of `input.xml`

