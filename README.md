---
layout: default
---
# Supplementary data for "Accurate nuclear quantum statistics on machine learned classical effective potentials"

Scripts, input files, codes and models required to reproduce the results of "Zaporozhets I. et al. Accurate nuclear quantum statistics on machine learned classical effective potentials
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

2) Download and install `te-pigs` with
```
pip install git+https://github.com/felixmusil/te-pigs.git
pip install git+https://github.com/felixmusil/mace.git@develop
```



4) Download and install `i-pi`: [https://github.com/venkatkapil24/i-pi.git](https://github.com/venkatkapil24/i-pi.git)

to use i-pi you will need to setup the environment variables using
    ```source $PATH2IPI/env.sh```

5) For bulk water only:
  Download  interface for MBPol calculations, `MBX-pythonic`:
  [https://github.com/venkatkapil24/MBX-pythonic.git](https://github.com/venkatkapil24/i-pi.git) and setup enfironmental variable `$MBX_HOME`, eg,
  `export $MBX_HOME=$HOME/software/MBX-pythonic/`



## License
The content of this repository is licensed under the CC-BY-SA-4.0 license. See the file `LICENSE` for details.