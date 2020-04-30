#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=96GB
#SBATCH --job-name=STARprefetch_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=prefetch_%j.out
  
module load edirect/20181128
module load sra-tools/intel/2.9.6
module load parallel/20171022
  
cd $SLURM_SUBMIT_DIR

parallel --verbose -j 20 prefetch --output-directory $SCRATCH/colobus_project/STAR/data/Rcolobus/ {} ::: $(grep -v Run SRR.numbers )
wait 

exit 

