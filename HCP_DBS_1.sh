#!/bin/bash

#Authored by Aaron Warren (aaron.warren@unimelb.edu.au)

#This code documents the image pre-processing steps relevant to the following paper:

#The optimal target and connectivity for DBS in Lennox-Gastaut syndrome
#Aaron E.L Warren, Linda J. Dalic, Kristian J. Bulluss, Annie Roten, Wesley Thevathasan, John S. Archer

#Analysis performed using the High Performance Computing (HPC) system ("Spartan") operated by Research Computing Services at The University of Melbourne:

#Lev Lafayette, Greg Sauter, Linh Vu, Bernard Meade, "Spartan Performance and Flexibility: An HPC-Cloud Chimera", OpenStack Summit, Barcelona, October 27, 2016. doi.org/10.4225/49/58ead90dceaaa
#https://dashboard.hpc.unimelb.edu.au 

#This bash script creates slurm submission scripts for each subject (*_sub.script). After per-subject scripts are created, they are ready to submit via sbatch.

#Usage: bash HCP_DBS_1.sh

#For each subject, the following steps are performed:

#1. Convert DWI data from .nii.gz to .mif and embed diffusion gradient encoding information within the image header
#2. Upsample DWI data to 1mm iso
#3. Bias correction of DWI data using ANTs
#4. Generate brain mask  
#5. Generate mean b-0 image
#6. Estimate subject-specific tissue response functions using 'dhollander' method (Dhollander, Thijs, David Raffelt, and Alan Connelly. "Unsupervised 3-tissue response function estimation from single-shell or multi-shell diffusion MR data without a co-registered T1 image." ISMRM Workshop on Breaking the Barriers of Diffusion MRI. Vol. 5. No. 5. ISMRM, 2016.)
#7. Copy the subject-specific tissue response functions to a common directory, for group averaging (see HCP_DBS_2.sh)

#load in a .txt file called "HCP_SUBJECTS.txt" listing HCP subject IDs (e.g., 100307, 100408, etc)

SUBJECTS=$(cat /data/HCP_SUBJECTS.txt);

#create slurm job for each subject 
for SUBJID in $SUBJECTS; do 
	
	#define a script name for each subject
	SCRIPTNAME=$SUBJID"_sub.script" 
	
	#define basic settings for slurm
	echo '#!/bin/bash' > $SCRIPTNAME
	echo '#SBATCH --time=01:0:0' >> $SCRIPTNAME
	echo '#SBATCH --ntasks=1' >> $SCRIPTNAME
	echo '#SBATCH --mem-per-cpu=32000' >> $SCRIPTNAME
	echo "#SBATCH --job-name=$SUBJID" >> $SCRIPTNAME
	echo >> $SCRIPTNAME 
		
	#load required software
	echo module load gcc/8.3.0 >> $SCRIPTNAME
	echo module load openmpi/3.1.4 >> $SCRIPTNAME
	echo module load mrtrix/3.0.1-python-2.7.16 >> $SCRIPTNAME

	#define/create some directories
	echo RAW_DIR=/data/HCP-RAW/${SUBJID}/T1w/Diffusion >> $SCRIPTNAME #where the minimally pre-processed ("raw") HCP data are stored
	echo PP_DIR=/data/HCP-PP/${SUBJID} >> $SCRIPTNAME #where the additionally processed data are stored
	echo RF_DIR=/data/HCP-PP/GROUP-RF >> $SCRIPTNAME #where the subject-specific tissue response functions should be copied to, for group averaging 
	echo "if [ ! -d ${PP_DIR} ]; then mkdir ${PP_DIR}; fi" >> $SCRIPTNAME
	echo "if [ ! -d ${RF_DIR} ]; then mkdir ${RF_DIR}; fi" >> $SCRIPTNAME
	
	#convert diffusion data to mrtrix's mif/non-compressed format and embed diffusion gradient encoding information within the image header
	echo mrconvert $RAW_DIR/data.nii.gz $PP_DIR/DWI.mif -fslgrad $RAW_DIR/bvecs $RAW_DIR/bvals >> $SCRIPTNAME 
	
	#upsample to 1mm
	echo mrgrid $PP_DIR/DWI.mif regrid $PP_DIR/DWI_1MM.mif -voxel 1 >> $SCRIPTNAME

	#bias field correction using ants
	echo dwibiascorrect ants $PP_DIR/DWI_1MM.mif $PP_DIR/DWI_1MM_BC.mif >> $SCRIPTNAME 
	
	#brain mask creation - this needs manual checking 
	echo dwi2mask $PP_DIR/DWI_1MM_BC.mif $PP_DIR/MASK.mif >> $SCRIPTNAME
	
	#create mean b-zero image and save as .nii.gz (used for various later steps - e.g., as reference image for warping MNI space ROIs to subject native DWI space)
	echo dwiextract $PP_DIR/DWI_1MM_BC.mif - -bzero | mrmath - mean $PP_DIR/MEANB0.nii.gz -axis 3 >> $SCRIPTNAME 
	
	#estimate subject-specific tissue response functions using 'dhollander' method
	echo dwi2response dhollander -voxels $PP_DIR/RF_VOXELS_DHOLL.mif $PP_DIR/DWI_1MM_BC.mif $PP_DIR/RF_WM_DHOLL.txt $PP_DIR/RF_GM_DHOLL.txt $PP_DIR/RF_CSF_DHOLL.txt
	
	#copy to group tissue response function directory
	echo cp $PP_DIR/RF_WM_DHOLL.txt $RF_DIR/$SUBJID_RF_WM_DHOLL.txt >> $SCRIPTNAME
	echo cp $PP_DIR/RF_GM_DHOLL.txt $RF_DIR/$SUBJID_RF_GM_DHOLL.txt >> $SCRIPTNAME
	echo cp $PP_DIR/RF_CSF_DHOLL.txt $RF_DIR/$SUBJID_RF_CSF_DHOLL.txt >> $SCRIPTNAME
	
done