#!/bin/bash
#SBATCH -J opt
#SBATCH --time=01:00:00
#SBATCH -n 1 -N1
#SBATCH --mem="20G"
#SBATCH --job-name=Morse
#SBATCH --partition=test

source ~/.bashrc
conda activate picg

python ../../../../scripts/linear_model_training.py --config config.yaml
