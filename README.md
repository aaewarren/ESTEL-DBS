# ESTEL-DBS

This code documents the image pre-processing and connectivity analysis steps relevant to the following paper:

The optimal target and connectivity for DBS in Lennox-Gastaut syndrome. Aaron E.L Warren, Linda J. Dalic, Kristian J. Bulluss, Annie Roten, Wesley Thevathasan, John S. Archer

Analyses were performed using the High Performance Computing (HPC) system ("Spartan") operated by Research Computing Services at The University of Melbourne:

https://dashboard.hpc.unimelb.edu.au
Lev Lafayette, Greg Sauter, Linh Vu, Bernard Meade, "Spartan Performance and Flexibility: An HPC-Cloud Chimera", OpenStack Summit, Barcelona, October 27, 2016. doi.org/10.4225/49/58ead90dceaaa

Analyses occur in two streams, one for the Human Connectome Project (HCP) DWI dataset (see scripts beginning HCP*) and another for the disease-matched LGS DWI dataset (see scripts beginning LGS*). The HCP and LGS streams differ slightly due to differences in DWI acquistion. In particular, for the LGS dataset, the following differences apply: (i) EPI distortion correction is performed without field-map or reverse phase-encoding data, using the ‘synthetic b0’ approach implemented by the excellent Synb0-DISCO software (https://github.com/MASILab/Synb0-DISCO); and (ii) FOD images are generated using Single-Shell 3-Tissue CSD (https://3tissue.github.io), rather than multi-shell multi-tissue CSD (used for the multi-shell HCP dataset). 

**HCP_DBS_1.sh and LGS_DBS_1.sh**

These bash files create a separate sub-script for each subject (*_sub.script). After per-subject scripts are created, they are ready to submit to SLURM via sbatch. They perform initial DWI and T1 pre-processing steps and calculate subject-specific tissue response functions using the 'dhollander' method in MRTrix3. 

**HCP_DBS_2.sh and LGS_DBS_2.sh**

These bash files are read to submit to SLURM. They perform group-averaging of subject-specific tissue response functions created above. 

**HCP_DBS_3.sh and LGS_DBS_3.sh**

These bash files create a separate sub-script for each subject (*_sub-2.script). After per-subject scripts are created, they are ready to submit to SLURM via sbatch. They create FOD images (using the group-average tissue response functions calculated in previous step) and perform tractography, SIFT2, and connectome calculation steps. Connectomes are saved as *csv files. 
