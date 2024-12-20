#!/bin/bash
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
module load R
Rscript 2E.2c-GO.R 

