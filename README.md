# ESTEL-DBS

This code documents the image pre-processing and connectivity analysis steps relevant to the following paper:

_The optimal target and connectivity for DBS in Lennox-Gastaut syndrome._
Aaron E.L Warren, Linda J. Dalic, Kristian J. Bulluss, Annie Roten, Wesley Thevathasan, John S. Archer

Analyses were performed using the High Performance Computing (HPC) system (["Spartan"](https://dashboard.hpc.unimelb.edu.au)) operated by Research Computing Services at The University of Melbourne:

Lev Lafayette, Greg Sauter, Linh Vu, Bernard Meade, "Spartan Performance and Flexibility: An HPC-Cloud Chimera", OpenStack Summit, Barcelona, October 27, 2016. doi.org/10.4225/49/58ead90dceaaa

Analyses occur in two streams, one for the Human Connectome Project (HCP) DWI dataset (see scripts beginning HCP*) and another for the disease-matched LGS DWI dataset (see scripts beginning LGS*). The HCP and LGS streams differ slightly due to differences in DWI acquistion. In particular, for the LGS dataset, the following differences apply: (i) EPI distortion correction is performed without field-map or reverse phase-encoding data, using the ‘synthetic b0’ approach implemented by the excellent [Synb0-DISCO software](https://github.com/MASILab/Synb0-DISCO); and (ii) FOD images are generated using [Single-Shell 3-Tissue CSD](https://3tissue.github.io), rather than multi-shell multi-tissue CSD (used for the multi-shell HCP dataset). 

_**HCP_DBS_1.sh and LGS_DBS_1.sh**_

These bash files create a separate sub-script for each subject (*_sub.script). After per-subject scripts are created, they are ready to submit to SLURM via sbatch. They perform initial DWI and T1 pre-processing steps and calculate subject-specific tissue response functions using the 'dhollander' method in MRTrix3. 

_**HCP_DBS_2.sh and LGS_DBS_2.sh**_

These bash files are read to submit to SLURM. They perform group-averaging of subject-specific tissue response functions created above. 

_**HCP_DBS_3.sh and LGS_DBS_3.sh**_

These bash files create a separate sub-script for each subject (*_sub-2.script). After per-subject scripts are created, they are ready to submit to SLURM via sbatch. They create FOD images (using the group-average tissue response functions calculated in previous step) and perform tractography, SIFT2, and connectome calculation steps. Connectomes are saved as *csv files. 
