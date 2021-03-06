#!/bin/bash

#Authored by Aaron Warren (aaron.warren@unimelb.edu.au)

#This code documents the image pre-processing steps relevant to the following paper:

#The optimal target and connectivity for DBS in Lennox-Gastaut syndrome
#Aaron E.L Warren, Linda J. Dalic, Kristian J. Bulluss, Annie Roten, Wesley Thevathasan, John S. Archer

#Analysis performed using the High Performance Computing (HPC) system ("Spartan") operated by Research Computing Services at The University of Melbourne:

#Lev Lafayette, Greg Sauter, Linh Vu, Bernard Meade, "Spartan Performance and Flexibility: An HPC-Cloud Chimera", OpenStack Summit, Barcelona, October 27, 2016. doi.org/10.4225/49/58ead90dceaaa
#https://dashboard.hpc.unimelb.edu.au 

#This bash script creates slurm submission scripts for each subject (*_sub-2.script). After per-subject scripts are created, they are ready to submit via sbatch.

#Usage: bash LGS_DBS_3.sh

#For each subject, the following steps are performed:

#1. Using T1 Freesurfer recon-all output calculated in LGS_DBS_1.sh, generate a "5 tissue type" segmentation using hybrid surface-volume method (Smith, R, Skoch, A, Bajada, CJ, et al. Hybrid surface-volume segmentation for improved anatomically-constrained tractography. Organisation for Human Brain Mapping. 2020:1034.)
#2. Align "5 tissue type" segmentation (computed above in step 1) to subject's T1 image, using the inverse of the DWI-to-T1 transform calculated in LGS_DBS_1.sh
#3. Generate fibre orientation distribution (FOD) images using single-shell 3-tissue CSD (https://3tissue.github.io/doc/ss3t-csd.html)
#4. 3-tissue bias field and intensity normalisation
#5. Use FSL's FIRST to create bilateral thalamic mask, for use in tracking, and then register this to DWI space (using the inverse of the DWI-to-T1 transform calculated in LGS_DBS_1.sh)
#6. Generate streamlines using both whole brain (10 million) and bilateral thalamic mask created in step 4 (10 million)
#7. Concatenate the whole-brain and bilateral thalamus tract files
#8. Perform SIFT2 on the concatenated tract file
#9. Generate connectomes (.csv) using whole-brain atlas which includes the per-ESTEL-patient bilateral VTAs (nonlinearly warped to each disease-matched LGS subject's native brain space)

#load in a .txt file called "LGS-SUBJECTS.txt" listing LGS subject IDs (e.g., 001, 002, etc)

SUBJECTS=$(cat /data/LGS_SUBJECTS.txt);

#create slurm job for each subject 
for SUBJID in $SUBJECTS; do 
	
	#define a script name for each subject
	SCRIPTNAME=$SUBJID"_sub-2.script" 
	
	#define basic settings for slurm
	echo '#!/bin/bash' > $SCRIPTNAME
	echo '#SBATCH --time=120:0:0' >> $SCRIPTNAME #note length of time - largely needed for streamline generation steps
	echo '#SBATCH --ntasks=1' >> $SCRIPTNAME
	echo '#SBATCH --mem-per-cpu=80000' >> $SCRIPTNAME #note memory usage - needed to combine the (very large) *.tck outputs
	echo "#SBATCH --job-name=$SUBJID" >> $SCRIPTNAME
	echo >> $SCRIPTNAME 
		
	#load required software
	echo module load gcc/8.3.0 >> $SCRIPTNAME
	echo module load openmpi/3.1.4 >> $SCRIPTNAME
	echo module load foss/2019b >> $SCRIPTNAME
	echo module load mrtrix3tissue/5.2.8-python-2.7.16 >> $SCRIPTNAME
	echo module load fsl/6.0.3-python-3.7.4 >> $SCRIPTNAME
	echo module load freesurfer/7.1.1-centos7_x86_64 >> $SCRIPTNAME
	
	#create/define some directories
	echo FS_DIR=/data/LGS-PP/${SUBJID}/FREESURFER7.1.1/${SUBJID} >> $SCRIPTNAME #where Freesurfer's recon-all output will be stored
	echo PP_DIR=/data/LGS-PP/${SUBJID} >> $SCRIPTNAME #where the additionally processed data are stored
	echo RF_DIR=/data/LGS-PP/GROUP-RF >> $SCRIPTNAME #where the group-averaged tissue response functions are stored
	echo VTA_DIR=/data/LGS-PP/${SUBJID}/VTAS >> $SCRIPTNAME #where the brain atlas files (ESTEL patient-specific bilateral VTAs), for use in connectome generation, are stored
	echo CONN_DIR=/data/LGS-PP/VTA_CONNECTOMES/${SUBJID} >> $SCRIPTNAME #where the connectome (*.csv) files should be copied, for subsequent analysis in MatLab (see LGS_DBS_4.m)
	echo FAST_DIR=/data/LGS-PP/${SUBJID}/FAST >> $SCRIPTNAME #where fsl FAST output/thalamic mask creation output is stored
	echo "if [ ! -d ${FAST_DIR} ]; then mkdir ${FAST_DIR}; fi" >> $SCRIPTNAME
	echo "if [ ! -d ${CONN_DIR} ]; then mkdir ${CONN_DIR}; fi" >> $SCRIPTNAME
	
	#generate a "5 tissue type" segmentation from freesurfer output using hybrid surface-volume method
	echo 5ttgen -nocrop hsvs $FS_DIR $PP_DIR/5tt_hsvs.mif >> $SCRIPTNAME
	
	#register "5 tissue type" segmentation to DWI - this step not performed for the HCP data because T1/Freesurfer data already aligned with DWI
	echo mrtransform $PP_DIR/5tt_hsvs.mif -linear $PP_DIR/DWI2T1_FSL2MRTRIX.txt -inverse -interp nearest $PP_DIR/5tt_hsvs2DWI.mif >> $SCRIPTNAME
	
	#generate fod images using single-shell 3-tissue CSD (https://3tissue.github.io/doc/ss3t-csd.html), making use of the group-averaged 3-tissue 'Dhollander'-estimated response functions
	echo ss3t_csd_beta1 $PP_DIR/DWI_DN_UR_EC_1MM_BC.mif $RF_DIR/GROUPAVG_RF_WM_DHOLL.txt $PP_DIR/WM_FOD_GROUP_RF.mif $RF_DIR/GROUPAVG_RF_GM_DHOLL.txt $PP_DIR/GM_FOD_GROUP_RF.mif $RF_DIR/GROUPAVG_RF_CSF_DHOLL.txt $PP_DIR/CSF_FOD_GROUP_RF.mif -mask $PP_DIR/MASK.mif >> $SCRIPTNAME

	#3-tissue bias field and intensity normalisation
	echo mtnormalise $PP_DIR/WM_FOD_GROUP_RF.mif $PP_DIR/WM_FOD_GROUP_RF_NORM.mif $PP_DIR/GM_FOD_GROUP_RF.mif $PP_DIR/GM_FOD_GROUP_RF_NORM.mif $PP_DIR/CSF_FOD_GROUP_RF.mif $PP_DIR/CSF_FOD_GROUP_RF_NORM.mif -mask $PP_DIR/MASK.mif >> $SCRIPTNAME

	#use FSL and MRTRix3 functions to create bilateral thalamus mask for use in tracking
	echo "cd ${FAST_DIR}" >> $SCRIPTNAME 
	echo mrconvert $FS_DIR/mri/norm.mgz T1.nii -stride -1,+2,+3 >> $SCRIPTNAME
	echo mrconvert $FS_DIR/mri/aparc+aseg.mgz aparc.mif >> $SCRIPTNAME
	echo run_first_all -s R_Thal,L_Thal -i T1.nii -b -o first >> $SCRIPTNAME
	echo meshconvert first-R_Thal_first.vtk first-R_Thal_transformed.vtk -transform first2real T1.nii >> $SCRIPTNAME
	echo mesh2voxel first-R_Thal_transformed.vtk aparc.mif Right-Thalamus-Proper.mif >> $SCRIPTNAME
	echo meshconvert first-L_Thal_first.vtk first-L_Thal_transformed.vtk -transform first2real T1.nii >> $SCRIPTNAME
	echo mesh2voxel first-L_Thal_transformed.vtk aparc.mif Left-Thalamus-Proper.mif >> $SCRIPTNAME
	echo mrtransform Right-Thalamus-Proper.mif -linear $PP_DIR/DWI2T1_FSL2MRTRIX.txt -inverse Right-Thalamus-Proper2DWI.nii.gz >> $SCRIPTNAME
	echo mrtransform Left-Thalamus-Proper.mif -linear $PP_DIR/DWI2T1_FSL2MRTRIX.txt -inverse Left-Thalamus-Proper2DWI.nii.gz >> $SCRIPTNAME
	echo fslmaths Right-Thalamus-Proper2DWI.nii.gz -thr 0.5 -bin $PP_DIR/Bin_Thr50_Right-Thalamus-Proper2DWI.nii.gz >> $SCRIPTNAME
	echo fslmaths Left-Thalamus-Proper2DWI.nii.gz -thr 0.5 -bin $PP_DIR/Bin_Thr50_Left-Thalamus-Proper2DWI.nii.gz >> $SCRIPTNAME
	echo fslmaths $PP_DIR/Bin_Thr50_Right-Thalamus-Proper2DWI.nii.gz -add $PP_DIR/Bin_Thr50_Left-Thalamus-Proper2DWI.nii.gz -bin $PP_DIR/Bin_Thr50_Bilat-Thalamus-Proper2DWI.nii.gz >> $SCRIPTNAME
	echo "cd $PP_DIR" >> $SCRIPTNAME 
	
	#generate whole-brain tractogram (10mill tracts) using iFOD2 algorithm (default in tckgen) and dynamic seeding - minimum lengh=2.5 mm and maximum length=250mm
	echo tckgen -select 10000000 -maxlength 250 -minlength 2.5 -act $PP_DIR/5tt_hsvs2DWI.mif -seed_dynamic $PP_DIR/WM_FOD_GROUP_RF_NORM.mif $PP_DIR/WM_FOD_GROUP_RF_NORM.mif $PP_DIR/tracts_10mill_wholebrain_ds_max250mm_act.tck >> $SCRIPTNAME
	
	#to ensure adequate inclusion of streamlines involving thalamic VTAs, generate additional tracts (10mill) using the bilateral thalamic mask created above
	echo tckgen -select 10000000 -maxlength 250 -minlength 2.5 -act $PP_DIR/5tt_hsvs2DWI.mif -seed_image $PP_DIR/Bin_Thr50_Bilat-Thalamus-Proper2DWI.nii.gz $PP_DIR/WM_FOD_GROUP_RF_NORM.mif $PP_DIR/tracts_10mill_bithal_seedimage_max250mm_act.tck >> $SCRIPTNAME
	
	#concatenate tract files into a single, 20mill tractogram 
	echo tckedit $PP_DIR/tracts_10mill_wholebrain_ds_max250mm_act.tck $PP_DIR/tracts_10mill_bithal_seedimage_max250mm_act.tck $PP_DIR/tracts_20mill_concat.tck >> $SCRIPTNAME
	
	#to save space, clean up large tract files we no longer need
	echo rm -f $PP_DIR/tracts_10mill_wholebrain_ds_max250mm_act.tck $PP_DIR/tracts_10mill_bithal_seedimage_max250mm_act.tck >> $SCRIPTNAME
	
	#perform sift2 on the concatenated file, saving mu parameter (subject-specific proportionality coefficient)
	echo tcksift2 -act $PP_DIR/5tt_hsvs2DWI.mif -out_mu $PP_DIR/sift2_mu.txt $PP_DIR/tracts_20mill_concat.tck $PP_DIR/WM_FOD_GROUP_RF_NORM.mif $PP_DIR/sift2_weights.txt >> $SCRIPTNAME
	
	#generate connectome (*.csv file) for each ESTEL patient, using whole-brain atlas which inlcludes the patient-specific bilateral VTA	
	echo "for vta in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; do
		
	tck2connectome $PP_DIR/tracts_20mill_concat.tck $VTA_DIR/${SUBJID}_ATLAS_429HCPMMPTIANCERES_VTA_${vta}.nii.gz $PP_DIR/connectome_sift2_vta_${vta}.csv -tck_weights_in $PP_DIR/sift2_weights.txt -keep_unassigned -assignment_radial_search 3;
		
	done" >> $SCRIPTNAME
	
	#clean up large concatenated tract file we no longer need
	echo rm -f $PP_DIR/tracts_20mill_concat.tck >> $SCRIPTNAME
	
	#copy the connectome *.csv files and sift2 mu output file to a group directory, for subsequent analysis in MatLab (see LGS_DBS_4.m)
	echo cp $PP_DIR/connectome_sift2_vta_${vta}.csv ${CONN_DIR} >> $SCRIPTNAME
	echo cp $PP_DIR/sift2_mu.txt ${CONN_DIR} >> $SCRIPTNAME
	
done