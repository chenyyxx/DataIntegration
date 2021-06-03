#!/bin/bash
#SBATCH -p shared
#SBATCH -t 0-20:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=32000



# load R
module load gcc/7.1.0-fasrc01 bedops/2.4.25-fasrc01
module load ucsc/20150820-fasrc01
module load Anaconda/5.0.1-fasrc02

conda create --name research python=3.8
conda activate research
conda install numpy
conda install pandas
pip install pyliftover
pip install pyBigWig
pip install dataintegrator