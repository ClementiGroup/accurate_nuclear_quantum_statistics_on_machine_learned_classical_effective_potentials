#!/bin/bash
#SBATCH -J opt
#SBATCH --time=2-00:00:00
#SBATCH -n 1 -N1
#SBATCH --mem="500M"
#SBATCH --job-name=Morse


source ~/.bashrc
conda activate picg

rm /tmp/ipi_morse_PIMD_efficient_sampling_100

if [ -f "RESTART" ]; then
    i-pi RESTART > log.i-pi &
else
	    i-pi input.xml > log.i-pi &
	    fi


sleep 30

for x in {1..64};
do
i-pi-driver -u -h morse_PIMD_efficient_sampling_100 -m morse &
done

wait

