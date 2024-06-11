#!/bin/bash
#SBATCH -J opt
#SBATCH --time=02:00:00
#SBATCH -n 1 -N1
#SBATCH --mem="5G"
#SBATCH --job-name=Morse
#SBATCH --partition=test

source ~/.bashrc
conda activate picg


python ../../../../scripts/generate_rescaled_dataset.py --config config.yaml
