#!/bin/sh

PBS_O_WORKDIR= /home/$USER/

cd $PBS_O_WORKDIR

module purge
module load gcc-4.6.2
module load mvapich2-1.9a2/gnu-4.6.2

echo "=== PART 1 ===" >> result.txt

mpicc -o part1 part1.c -lm

P1N1=$(qsub JobFile_P1_N1)
P1N2=$(qsub -W depend=afterany:${P1N1} JobFile_P1_N2)
P1N4=$(qsub -W depend=afterany:${P1N2} JobFile_P1_N4)
qsub -W depend=afterany:${P1N4} JobFile_P1_N8
