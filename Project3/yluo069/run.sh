#!/bin/sh

PBS_O_WORKDIR= /home/$USER/

cd $PBS_O_WORKDIR

module purge
module load gcc-4.6.2
module load mvapich2-1.9a2/gnu-4.6.2

echo "=== Project 3 ===" >> results

mpicc -o part1 part1.c -lm
mpicc -o part2 part2.c -lm
mpicc -o part3 part3.c -lm

P1N1=$(qsub JobFile_P1_N1)
P1N2=$(qsub -W depend=afterany:${P1N1} JobFile_P1_N2)
P1N4=$(qsub -W depend=afterany:${P1N2} JobFile_P1_N4)
P1N8=$(qsub -W depend=afterany:${P1N4} JobFile_P1_N8)

P2N1=$(qsub -W depend=afterany:${P1N8} JobFile_P2_N1)
P2N2=$(qsub -W depend=afterany:${P2N1} JobFile_P2_N2)
P2N4=$(qsub -W depend=afterany:${P2N2} JobFile_P2_N4)
P2N8=$(qsub -W depend=afterany:${P2N4} JobFile_P2_N8)
P3N1=$(qsub -W depend=afterany:${P2N8} JobFile_P3_N1)
P3N2=$(qsub -W depend=afterany:${P3N1} JobFile_P3_N2)
P3N4=$(qsub -W depend=afterany:${P3N2} JobFile_P3_N4)

qsub -W depend=afterany:${P3N4} JobFile_P3_N8
