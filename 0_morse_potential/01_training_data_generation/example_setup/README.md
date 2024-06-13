---
layout: default
---

# PIMD simulation

Template files for PIMD simulation. 

Things to change in input files for running the simultation:
input.xml:  line  20   : Set correct temperature instead of 100 
            line  16   : Use temperature as a suffix for socket address (replace 100 with correct temperature)

submit.sh: line 7      : Set correct suffix for the file to be removed  (replace 100 with correct temperature)
           line 15     : Set the same socket address as in line 16 of input.xml
