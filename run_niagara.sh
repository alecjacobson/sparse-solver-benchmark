#!/bin/sh

#SBATCH --cpus-per-task=40
#SBATCH --export=ALL
#SBATCH --job-name="sym"
#SBATCH --nodes=1
#SBATCH --output="fusy.%j.%N.out"
#SBATCH -t 2:00:00

module load NiaEnv/2019b
module load cmake/3.17.3
#module load intel
module load intel/2019u4
#module load intel
#module load gcc
module load metis/5.1.0
source ${MKLROOT}/bin/mklvars.sh intel64


export OMP_NUM_THREADS=20
export MKL_NUM_THREADS=20



cmake ../ -DCMAKE_PREFIX_PATH=$MKLROOT/lib/intel64  -DCMAKE_BUILD_TYPE=Release
cmake ../ -DCMAKE_PREFIX_PATH=$MKLROOT/include  -DCMAKE_BUILD_TYPE=Release
make

./sparse_solver_benchmark ../xyzrgb_dragon-720K.ply



