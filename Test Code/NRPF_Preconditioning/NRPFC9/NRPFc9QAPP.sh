#!/bin/bash
#SBATCH --job-name=NRPF_case9QAPP    # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dj.neufeld@uleth.ca     # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
pwd; hostname; date

module load python/3.9
module load scipy-stack
module load cmake
module load symengine
source /home/dneufeld/env/bin/activate

echo "Running..."

python NRPFCase9QAPP.py

date

deactivate
