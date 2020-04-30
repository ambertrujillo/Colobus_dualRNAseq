#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=24GB
#SBATCH --job-name=STARdump_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aet359@nyu.edu
#SBATCH --output=dump_%j.out
  
module load edirect/20181128
module load sra-tools/intel/2.9.6
module load parallel/20171022
  
cd $SLURM_SUBMIT_DIR

parallel --verbose -j 20 fastq-dump --split-files --outdir $SCRATCH/colobus_project/STAR/data/Rcolobus/ {} ::: $(ls $SCRATCH/colobus_project/STAR/data/Rcolobus/*.sra )
wait 

exit 
