#!/bin/bash

#Authored by Aaron Warren (aaron.warren@unimelb.edu.au)

#This code documents the image pre-processing steps relevant to the following paper:

#The optimal target and connectivity for DBS in Lennox-Gastaut syndrome
#Aaron E.L Warren, Linda J. Dalic, Kristian J. Bulluss, Annie Roten, Wesley Thevathasan, John S. Archer

#Analysis performed using the High Performance Computing (HPC) system ("Spartan") operated by Research Computing Services at The University of Melbourne:

#Lev Lafayette, Greg Sauter, Linh Vu, Bernard Meade, "Spartan Performance and Flexibility: An HPC-Cloud Chimera", OpenStack Summit, Barcelona, October 27, 2016. doi.org/10.4225/49/58ead90dceaaa
#https://dashboard.hpc.unimelb.edu.au 

#This bash script performs group-level averaging of subject-specific tissue response functions and is ready to submit via sbatch.

#Usage: sbatch LGS_DBS_2.sh

#SBATCH --time=01:0:0
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=LGS_DBS_2

#load required software
module load gcc/8.3.0
module load openmpi/3.1.4
module load mrtrix/3.0.1-python-2.7.16

#cd to directory and perform group averaging of subject-specific tissue response functions

cd /data/LGS-PP/GROUP-RF
responsemean *_WM_DHOLL.txt GROUPAVG_RF_WM_DHOLL.txt
responsemean *_GM_DHOLL.txt GROUPAVG_RF_GM_DHOLL.txt
responsemean *_CSF_DHOLL.txt GROUPAVG_RF_CSF_DHOLL.txt