---
layout: default
---

# PIMD simulation

Template files to run one server for each simulation. 

Parameters to change in input files for running the simultation:
`input.xml`:  line  21   : Set correct temperature instead of 100 
            line  17   : Use temperature as a suffix for socket address (replace 100 with correct temperature)

`submit.sh`: lines 14-15      : Set correct suffix for the file to be removed  (replace 100 with correct temperature)

`run-driver.sh`:
           line 1   : Set the same socket address as in line 17 of `input.xml`


