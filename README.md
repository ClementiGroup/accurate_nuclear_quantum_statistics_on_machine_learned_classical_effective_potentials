---
layout: default
---
# Supplementary data for "Accurate nuclear quantum statistics at the cost of classical molecular dynamics"

Scripts, input files, codes and models required to reproduce the results of "Zaporozhets I. et al. Accurate nuclear quantum statistics at the cost of classical molecular
dynamics (2024)"

## Systems


- `0_morse_potential` : particle in 1D Morse potential
- `1_h2o_molecule` : single H2O molecule in vacuum
- `2_zundel_cation` : Zundel cation
- `3_bulk_h2o`:  Water box (256 molecules)


## Setup instructions

1) Setup conda environment:

`cd setup`
`./setup_environment.sh`

2) Download and install `te-pigs`


3) Setup an environment variable, pointing to the location of `te-pigs` package, for example:

`export TEPIGS_PATH=/home/username/software/te-pigs/`

4) Download and install `i-pi`: `https://github.com/venkatkapil24/i-pi.git`


## License
The content of this repository is licensed under the CC-BY-SA-4.0 license. See the file `LICENSE` for details.