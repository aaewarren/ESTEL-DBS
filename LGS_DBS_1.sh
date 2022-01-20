#!/bin/bash

#Authored by Aaron Warren (aaron.warren@unimelb.edu.au)

#This code documents the image pre-processing steps relevant to the following paper:

#The optimal target and connectivity for DBS in Lennox-Gastaut syndrome
#Aaron E.L Warren, Linda J. Dalic, Kristian J. Bulluss, Wesley Thevathasan, Annie Roten, John S. Archer

#Analysis performed using the High Performance Computing (HPC) system ("Spartan") operated by Research Computing Services at The University of Melbourne:

#Lev Lafayette, Greg Sauter, Linh Vu, Bernard Meade, "Spartan Performance and Flexibility: An HPC-Cloud Chimera", OpenStack Summit, Barcelona, October 27, 2016. doi.org/10.4225/49/58ead90dceaaa
#https://dashboard.hpc.unimelb.edu.au 

#This bash script creates slurm submission scripts for each subject (*_sub.script). After per-subject scripts are created, they are ready to submit via sbatch.

#Usage: bash LGS_DBS_1.sh

#For each subject, the following steps are performed:

#1. Perform initial DWI pre-processing steps (denoising, unringing, mask creation)
#2. Prepare inputs for Synb0-DISCO (see https://github.com/MASILab/Synb0-DISCO), and then run via docker 
#3. Submit outputs of Synb0-DISCO to FSL's eddy (see https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy)
#4. Upsample DWI to 1mm iso
#5. Bias correction of DWI data using ANTs
#6. Generate brain mask  
#7. Generate mean b0 image
#8. Estimate subject-specific tissue response functions using 'dhollander' method (Dhollander, Thijs, David Raffelt, and Alan Connelly. "Unsupervised 3-tissue response function estimation from single-shell or multi-shell diffusion MR data without a co-registered T1 image." ISMRM Workshop on Breaking the Barriers of Diffusion MRI. Vol. 5. No. 5. ISMRM, 2016.)
#9. Copy the subject-specific tissue response functions to a common directory, for group averaging (see LGS_DBS_2.sh)
#10. Run Freesurfer's recon-all (version 7.1.1) on the T1 image
#11. Create head- and brain-masked T1 image (using Freesurfer output) for use in T1-to-DWI co-registration
#12. Calculate DWI-to-T1 transform using FSL's epi_reg (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/UserGuide#epi_reg)
#13. Convert FSL-style DWI-to-T1 transform to MRTrix-style transform 
#14. Co-register T1 to DWI using inverse of the DWI-to-T1 transform 

#load in a .txt file called "LGS-SUBJECTS.txt" listing LGS subject IDs (e.g., 001, 002, etc)

SUBJECTS=$(cat /data/LGS-SUBJECTS.txt);

#create slurm job for each subject 
for SUBJID in $SUBJECTS; do 
	
	#define a script name for each subject
	SCRIPTNAME=$SUBJID"_sub.script" 
	
	#define basic settings for slurm
	echo '#!/bin/bash' > $SCRIPTNAME
	echo '#SBATCH --time=50:0:0' >> $SCRIPTNAME #note extra time for SynB0-DISCO and Freesurfer recon-all steps 
	echo '#SBATCH --ntasks=1' >> $SCRIPTNAME
	echo '#SBATCH --mem-per-cpu=32000' >> $SCRIPTNAME
	echo "#SBATCH --job-name=$SUBJID" >> $SCRIPTNAME
	echo >> $SCRIPTNAME 
		
	#load required software
	echo module load gcc/8.3.0 >> $SCRIPTNAME
	echo module load openmpi/3.1.4 >> $SCRIPTNAME
	echo module load mrtrix/3.0.1-python-2.7.16 >> $SCRIPTNAME
	echo module load fsl/6.0.1-python-3.7.4 >> $SCRIPTNAME
	echo module load foss/2019b >> $SCRIPTNAME #required for ANTs
	echo module load ants/2.3.2-python-3.7.4 >> $SCRIPTNAME
	echo module load freesurfer/7.1.1-centos7_x86_64 >> $SCRIPTNAME
	echo "source /usr/local/easybuild-2019/easybuild/software/core/freesurfer/7.1.1-centos7_x86_64/SetUpFreeSurfer.sh" >> $SCRIPTNAME #source freesurfer setup file
	
	#define/create some directories
	echo RAW_DIR=/data/LGS-RAW/${SUBJID} >> $SCRIPTNAME #where the raw DWI (converted from DICOM to *.mif format) and T1 data (in *.nii.gz format) are stored
	echo PP_DIR=/data/LGS-PP/${SUBJID} >> $SCRIPTNAME #where the pre-processed data are stored
	echo SYNB0_DIR=/data/LGS-PP/${SUBJID}/Synb0-DISCO >> $SCRIPTNAME #where the Synb0-DISCO inputs/outputs are stored
	echo SYNB0_DIR_IN=/data/LGS-PP/${SUBJID}/Synb0-DISCO/INPUTS >> $SCRIPTNAME #where the Synb0-DISCO inputs are stored
	echo SYNB0_DIR_OUT=/data/LGS-PP/${SUBJID}/Synb0-DISCO/OUTPUTS >> $SCRIPTNAME #where the Synb0-DISCO outputs are stored
	echo FS_DIR=/data/LGS-PP/${SUBJID}/FREESURFER7.1.1 >> $SCRIPTNAME #where Freesurfer's recon-all output will be stored
	echo RF_DIR=/data/LGS-PP/GROUP-RF >> $SCRIPTNAME #where the subject-specific tissue response functions should be copied to, for group averaging
	echo "if [ ! -d ${PP_DIR} ]; then mkdir ${PP_DIR}; fi" >> $SCRIPTNAME
	echo "if [ ! -d ${SYNB0_DIR} ]; then mkdir ${SYNB0_DIR}; fi" >> $SCRIPTNAME
	echo "if [ ! -d ${SYNB0_DIR_IN} ]; then mkdir ${SYNB0_DIR_IN}; fi" >> $SCRIPTNAME
	echo "if [ ! -d ${SYNB0_DIR_OUT} ]; then mkdir ${SYNB0_DIR_OUT}; fi" >> $SCRIPTNAME
	echo "if [ ! -d ${FS_DIR} ]; then mkdir ${FS_DIR}; fi" >> $SCRIPTNAME
	echo "if [ ! -d ${RF_DIR} ]; then mkdir ${RF_DIR}; fi" >> $SCRIPTNAME
	
	#perform initial denoising, unringing, and DWI mask creation (from distorted images) steps 
	echo dwidenoise ${RAW_DIR}/DWI.mif ${PP_DIR}/DWI_DN.mif >> $SCRIPTNAME
	echo mrdegibbs ${PP_DIR}/DWI_DN.mif ${PP_DIR}/DWI_DN_UR.nii.gz >> $SCRIPTNAME #note save as .nii.gz - as this will be used for Synb0-DISCO below
	echo dwi2mask ${RAW_DIR}/DWI.mif ${PP_DIR}/DWI_MASK_DIST.nii.gz >> $SCRIPTNAME #note save as .nii.gz - as this will be used for FSL's eddy command below
	
	#calculate/copy inputs required for Synb0-DISCO
	echo mrinfo ${RAW_DIR}/DWI.mif -export_grad_fsl ${SYNB0_DIR_OUT}/bvecs.txt ${SYNB0_DIR_OUT}/bvals.txt >> $SCRIPTNAME #bvals and bvecs 
	echo dwiextract ${RAW_DIR}/DWI.mif - -bzero | mrmath - mean ${SYNB0_DIR_IN}/b0.nii.gz -axis 3 >> $SCRIPTNAME #mean b0 image as .nii.gz 
	echo cp ${PP_DIR}/DWI_DN_UR.nii.gz ${SYNB0_DIR_IN}/DWI.nii.gz >> $SCRIPTNAME #denoised, unringed DWI images
	echo cp ${RAW_DIR}/T1.nii.gz ${SYNB0_DIR_IN} >> $SCRIPTNAME #T1 image
	echo cp /data/LGS-PP/license.txt ${SYNB0_DIR} >> #freesurfer licence file 
	
	#note also that a file called "acqparams.txt" must also be saved in ${SYNBO_DIR_IN}. i suggest creating this file manually for each patient to ensure correct params.
	#from the SynB0-DISCO help page: 
	#acqparams.txt describes the acqusition parameters, and is described in detail on the FslWiki for topup (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup).Briefly,
	#it describes the direction of distortion and tells TOPUP that the synthesized image
	#has an effective echo spacing of 0 (infinite bandwidth). An example acqparams.txt is
	#displayed below, in which distortion is in the second dimension, note that the second
	#row corresponds to the synthesized, undistorted, b0:
	#$ cat acqparams.txt 
	#0 1 0 0.062
	#0 1 0 0.000
	
	#now run Synb0-DISCO via docker. note that a freesurfer license file (license.txt) must be in /data/LGS-PP 
	echo "cd ${SYNB0_DIR}" >> $SCRIPTNAME
	echo docker run --rm -v $(pwd)/INPUTS/:/INPUTS/ -v $(pwd)/OUTPUTS:/OUTPUTS/ -v /data/LGS-PP/license.txt:/extra/freesurfer/license.txt --user $(id -u):$(id -g) hansencb/synb0 >> $SCRIPTNAME
	
	#after Synb0-DISCO finishes, run FSL's eddy, then convert output (i.e., unwarped images) back to *.mif format with embedded bvecs/bvals, for subsequent processing in mrtrix3
	echo "cd ${SYNB0_DIR_OUT}" >> $SCRIPTNAME
	echo eddy_openmp --imain=${PP_DIR}/DWI_DN_UR.nii.gz \
		--mask=${PP_DIR}/DWI_MASK_DIST.nii.gz \
		--acqp=${SYNBO_DIR_IN}/acqparams.txt \
		--index=index.txt --bvecs=bvecs.txt --bvals=bvals.txt --topup=topup --out=eddy_unwarped_images -v >> $SCRIPTNAME 
	echo mrconvert eddy_unwarped_images.nii.gz -fslgrad eddy_unwarped_images.eddy_rotated_bvecs bvals.txt ${PP_DIR}/DWI_DN_UR_EC.mif >> $SCRIPTNAME #EC = eddy correction
	echo "cd ${PP_DIR}" >> $SCRIPTNAME
	
	#upsample to 1mm
	echo mrgrid $PP_DIR/DWI_DN_UR_EC.mif regrid $PP_DIR/DWI_DN_UR_EC_1MM.mif -voxel 1 >> $SCRIPTNAME
	
	#bias correction in ants
	echo dwibiascorrect ants $PP_DIR/DWI_DN_UR_EC_1MM.mif $PP_DIR/DWI_DN_UR_EC_1MM_BC.mif >> $SCRIPTNAME
	
	#brain mask creation (unwarped, 1mm images)- this needs manual checking 
	echo dwi2mask $PP_DIR/DWI_DN_UR_EC_1MM_BC.mif $PP_DIR/MASK.mif >> $SCRIPTNAME

	#create mean b-zero image and save as .nii.gz (used for various later steps - e.g., as reference image for warping MNI space ROIs to subject native DWI space)
	echo dwiextract $PP_DIR/DWI_DN_UR_EC_1MM_BC.mif - -bzero | mrmath - mean $PP_DIR/MEANB0.nii.gz -axis 3 >> $SCRIPTNAME
	
	#estimate subject-specific tissue response functions using 'dhollander' method
	echo dwi2response dhollander -voxels $PP_DIR/RF_VOXELS_DHOLL.mif $PP_DIR/DWI_DN_UR_EC_1MM_BC.mif $PP_DIR/RF_WM_DHOLL.txt $PP_DIR/RF_GM_DHOLL.txt $PP_DIR/RF_CSF_DHOLL.txt >> $SCRIPTNAME
	
	#copy to group tissue response function directory
	echo cp $PP_DIR/RF_WM_DHOLL.txt $RF_DIR/$SUBJID_RF_WM_DHOLL.txt >> $SCRIPTNAME
	echo cp $PP_DIR/RF_GM_DHOLL.txt $RF_DIR/$SUBJID_RF_GM_DHOLL.txt >> $SCRIPTNAME
	echo cp $PP_DIR/RF_CSF_DHOLL.txt $RF_DIR/$SUBJID_RF_CSF_DHOLL.txt >> $SCRIPTNAME
	
	#now run Freesurfer's recon-all on the subject's T1 image
	echo SUBJECTS_DIR=${FS_DIR} >> $SCRIPTNAME #set SUBJECTS_DIR
	echo recon-all -i ${RAW_DIR}/T1.nii.gz -subjid $SUBJID -all >> $SCRIPTNAME	
	
	#create head and brain masks of T1 image using Freesurfer output, and apply some bias correction in ANTs
	echo mrconvert ${SUBJECTS_DIR}/$SUBJID/mri/orig/001.mgz $PP_DIR/origT1.nii.gz -stride +1,+2,+3 -force >> $SCRIPTNAME	
	echo bet $PP_DIR/origT1.nii.gz $PP_DIR/bet -A >> $SCRIPTNAME #-A option to get skull/scalp surfaces - to create whole-head mask 
	echo rm -f $PP_DIR/bet*mesh* $PP_DIR/bet.nii.gz $PP_DIR/bet*skull* >> $SCRIPTNAME
	echo mv $PP_DIR/bet_outskin_mask.nii.gz $PP_DIR/T1_headmask.nii.gz >> $SCRIPTNAME
	echo fslmaths $PP_DIR/origT1.nii.gz -mul $PP_DIR/T1_headmask.nii.gz $PP_DIR/T1_head.nii.gz >> $SCRIPTNAME
	echo N4BiasFieldCorrection -d 3 -i $PP_DIR/T1_head.nii.gz -x $PP_DIR/T1_headmask.nii.gz -o $PP_DIR/T1_head_BC.nii.gz >> $SCRIPTNAME
	echo mri_label2vol --seg ${SUBJECTS_DIR}/$SUBJID/mri/aparc+aseg.mgz --temp ${SUBJECTS_DIR}/$SUBJID/mri/rawavg.mgz --o $PP_DIR/aparc+aseg2rawavg.mgz --regheader ${SUBJECTS_DIR}/$SUBJID/mri/aparc+aseg.mgz >> $SCRIPTNAME
	echo mrconvert $PP_DIR/aparc+aseg2rawavg.mgz $PP_DIR/aparc+aseg.nii.gz -stride +1,+2,+3 -force >> $SCRIPTNAME
	echo fslmaths $PP_DIR/aparc+aseg.nii.gz -bin -kernel sphere 3 -dilM -ero -dilM -ero -dilM -ero -binv $PP_DIR/T1_mask_tmp1.nii.gz >> $SCRIPTNAME
	echo connectedcomp $PP_DIR/T1_mask_tmp1.nii.gz $PP_DIR/T1_mask_tmp2.nii.gz >> $SCRIPTNAME
	echo fslmaths $PP_DIR/T1_mask_tmp2.nii.gz -uthr 1 -binv $PP_DIR/T1_brainmask.nii.gz >> $SCRIPTNAME
	echo rm -f $PP_DIR/T1_mask_tmp?.nii.gz $PP_DIR/aparc+aseg2rawavg.mgz >> $SCRIPTNAME
	echo fslmaths $PP_DIR/T1_head_BC.nii.gz -mul $PP_DIR/T1_brainmask.nii.gz $PP_DIR/T1_brain_BC.nii.gz >> $SCRIPTNAME
	
	#now perform DWI-to-T1 coregistration using FSL's epi_reg procedure (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/UserGuide#epi_reg), using whole-head and brain-extracted T1 images computed above
	echo "cd $PP_DIR" >> $SCRIPTNAME
	echo epi_reg --epi=MEANB0.nii.gz --t1=T1_head_BC.nii.gz --t1brain=T1_brain_BC.nii.gz --out=DWI2T1 >> $SCRIPTNAME
	
	#convert fsl-style transformation matrix to mrtrix-style matrix
	echo transformconvert $PP_DIR/DWI2T1.mat $PP_DIR/MEANB0.nii.gz $PP_DIR/T1_head_BC.nii.gz flirt_import $PP_DIR/DWI2T1_FSL2MRTRIX.txt >> $SCRIPTNAME
	
	#register T1 to DWI using inverse of transform calculated by epi_reg - note the "-inverse" flag
	echo mrtransform $PP_DIR/T1_head_BC.nii.gz -linear $PP_DIR/DWI2T1_FSL2MRTRIX.txt -inverse $PP_DIR/T1_head_BC_2DWI.nii.gz >> $SCRIPTNAME
		
done