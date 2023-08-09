#!/bin/bash
#SBATCH --job-name=NRPF_case3    # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dj.neufeld@uleth.ca     # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=8gb                     # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
pwd; hostname; date

module load python/3.9
module load scipy-stack
module load cmake
module load symengine
source /home/dneufeld/env/bin/activate

echo "Running..."

python NRPFCase3.py

date

deactivate
