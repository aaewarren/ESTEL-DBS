# ESTEL-DBS

This code documents the image pre-processing and connectivity analysis steps relevant to the following paper:

_The optimal target and connectivity for DBS in Lennox-Gastaut syndrome._
Aaron E.L Warren, Linda J. Dalic, Kristian J. Bulluss, Annie Roten, Wesley Thevathasan, John S. Archer

Analyses were performed using the High Performance Computing (HPC) system ([Spartan](https://dashboard.hpc.unimelb.edu.au)) operated by Research Computing Services at The University of Melbourne:

Lev Lafayette, Greg Sauter, Linh Vu, Bernard Meade, "Spartan Performance and Flexibility: An HPC-Cloud Chimera", OpenStack Summit, Barcelona, October 27, 2016. doi.org/10.4225/49/58ead90dceaaa

Analyses occur in two streams, one for the Human Connectome Project (HCP) DWI dataset (scripts beginning HCP*) and another for the disease-matched LGS DWI dataset (scripts beginning LGS*). The HCP and LGS streams differ slightly due to differences in DWI acquistion. In particular, for the LGS dataset, the following differences apply: (i) echo-planar imaging (EPI) distortion correction is performed without field-map or reverse phase-encoding data, using the ‘synthetic b0’ approach deployed by the excellent [Synb0-DISCO software](https://github.com/MASILab/Synb0-DISCO); and (ii) fibre orientation distribution (FOD) images are created using [Single-Shell 3-Tissue CSD](https://3tissue.github.io), rather than multi-shell multi-tissue CSD (as used for the multi-shell HCP dataset). 

Each stream includes four stages, executed consecutively (three bash scripts: *_1.sh, *_2.sh, *_3.sh; and lastly one MatLab script: *_4.m).

### **HCP_DBS_1.sh/LGS_DBS_1.sh**

These bash scripts create a separate sub-script for each subject (*_sub.script). After per-subject scripts are created, they are ready to submit to SLURM (via sbatch). They perform initial DWI and T1 pre-processing steps and calculate subject-specific tissue response functions using the [Dhollander](https://mrtrix.readthedocs.io/en/latest/reference/commands/dwi2response.html#dwi2response-dhollander) approach in [MRtrix3](https://mrtrix.readthedocs.io/en/latest/index.html) software. 

### **HCP_DBS_2.sh/LGS_DBS_2.sh**

These bash scripts are read to submit to SLURM (via sbatch). They perform [group-averaging](https://mrtrix.readthedocs.io/en/latest/reference/commands/responsemean.html) of the subject-specific tissue response functions created above. 

### **HCP_DBS_3.sh/LGS_DBS_3.sh**

These bash scripts create a separate sub-script for each subject (*_sub-2.script). After per-subject scripts are created, they are ready to submit to SLURM (via sbatch). They create FOD images (using the group-averaged tissue response functions calculated in previous step) and perform [tractography](https://mrtrix.readthedocs.io/en/latest/reference/commands/tckgen.html), [SIFT2](https://mrtrix.readthedocs.io/en/latest/reference/commands/tcksift2.html), and [connectome calculation](https://mrtrix.readthedocs.io/en/latest/reference/commands/tck2connectome.html) steps. Connectomes are saved as .csv files, for importing into MatLab below.

### **HCP_DBS_4.m/LGS_DBS_4.m**

These MatLab scripts perform post-processing steps on the structural connectome files (*.csv) computed by HCP_DBS_3.sh/LGS_DBS_3.sh, and then calculate parcel-wise correlations between connectivity values and ESTEL patient clinical outcomes. The scripts require the function "nii_tool.m" for reading/writing NIfTI data. This function can be downloaded from [here](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/42997/versions/80/previews/nii_tool.m/index.html?access_key=).

