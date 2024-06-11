#!/bin/bash
#SBATCH -J opt
#SBATCH --time=336:00:00
#SBATCH -n 16 -N1
#SBATCH --mem="500M"
#SBATCH --job-name=h5o2+_PIMD_efficient_sampling_0100K

source ~/.bashrc
conda activate picg

ipi=i-pi

export OMP_NUM_THREADS=1

touch /tmp/ipi_h5o2+_PIMD_efficient_sampling_0100K
rm /tmp/ipi_h5o2+_PIMD_efficient_sampling_0100K

# runs i-PI + driver to accumulate data

if [ -f "RESTART" ]; then
${ipi} RESTART > log.i-pi &
else
${ipi} input.xml > log.i-pi &
fi

sleep 30

for x in {1..16};
do
bash run-driver.sh &
done

wait
