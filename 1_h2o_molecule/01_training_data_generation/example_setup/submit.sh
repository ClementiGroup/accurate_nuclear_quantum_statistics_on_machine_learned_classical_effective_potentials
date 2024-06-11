#!/bin/bash
#SBATCH -J opt
#SBATCH --time=336:00:00
#SBATCH -n 1 -N1
#SBATCH --mem="100M"
#SBATCH --job-name=h2o-molecule_PIMD_efficient_sampling_0100K

source ~/.bashrc
conda activate picg

export OMP_NUM_THREADS=1


touch /tmp/ipi_h2o-molecule_PIMD_efficient_sampling_0100K
rm /tmp/ipi_h2o-molecule_PIMD_efficient_sampling_0100K

# runs i-PI + driver to accumulate data

if [ -f "RESTART" ]; then
i-pi RESTART > log.i-pi &
else
i-pi input.xml > log.i-pi &
fi

sleep 30

for x in {1..1};
do
bash run-driver.sh &
done

wait
