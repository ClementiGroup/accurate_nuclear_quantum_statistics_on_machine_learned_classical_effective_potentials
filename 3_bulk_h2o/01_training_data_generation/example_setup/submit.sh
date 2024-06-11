#!/bin/bash

#SBATCH --job-name=mbpol-test            # Job name, will show up in squeue output
#SBATCH --partition=main                # main or test
#SBATCH --ntasks=8       # Number of cores
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1                      # Ensure that all cores are on one machine
#SBATCH --time=5-00:30:00              # Runtime in DAYS-HH:MM:SS format
#SBATCH --mem="20GB"              # Memory per cpu in MB (see also --mem)
#SBATCH --output=job.out           # File to which standard out will be written
#SBATCH --error=job.err            # File to which standard err will be written


source ~/.bashrc
conda activate picg


ipi=i-pi
driver=$MBX_HOME/plugins/i-pi/bin/driver

rm /tmp/ipi_h2o_MBPOL-liquid

if [ -f "RESTART" ]; then
    i-pi RESTART > log.i-pi &
else
    i-pi input.xml > log.i-pi &
fi


sleep 30

export OMP_NUM_THREADS=2

for x in {1..8}
do
${driver} config.nrg mbx.json &
done

wait

