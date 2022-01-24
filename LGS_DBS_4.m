%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Authored by Aaron Warren (aaron.warren@unimelb.edu.au)

%This code documents the image processing steps relevant to the following paper:

%The optimal target and connectivity for DBS in Lennox-Gastaut syndrome
%Aaron E.L Warren, Linda J. Dalic, Kristian J. Bulluss, Annie Roten, Wesley Thevathasan, John S. Archer

%Analysis performed using the High Performance Computing (HPC) system ("Spartan") operated by Research Computing Services at The University of Melbourne:

%Lev Lafayette, Greg Sauter, Linh Vu, Bernard Meade, "Spartan Performance and Flexibility: An HPC-Cloud Chimera", OpenStack Summit, Barcelona, October 27, 2016. doi.org/10.4225/49/58ead90dceaaa
%https://dashboard.hpc.unimelb.edu.au 

%This MatLab script performs post-processing steps on structrual
%connectivity matrices output from MRtrix3 (*.csv), as well as correlations with clinical outcomes (Diaries, EEG, Diary-EEG Average)

%Usage: run LGS_DBS_4.m

%The following steps are performed:

%1. Load in per-LGS-subject/per-ESTEL-patient-VTA structural connectivity matrices (*.csv files), computed in LGS_DBS_3.sh
%2. Retain columns/rows of interest (i.e., left and right VTAs to all other brain regions within the 429-region parcellation)
%3. Multiply each matrix element by subject-specific fibre-density proportionality coefficient (mu), output from SIFT2 (see LGS_DBS_3.sh)
%4. Average across LGS subjects (yielding per-ESTEL-patient connectivity values, for left and right VTAs separately)
%5. Sum together LGS subject-averaged matrices corresponding to left and right VTAs (yielding per-ESTEL-patient bilateral connectivity values, for left and right VTAs combined)
%6. Average across the per-ESTEL-patient bilateral connectivity values (yielding a single 429-length vector representing average bilteral VTA connectivity across all ESTEL patients)
%7. Square-root normalise the connectivity vector from step 6 for display purposes, and save as .nii
%8. Run parcel-wise correlations (connectivity vs % seizure reduction), separately for each outcome measure (Diaries, EEG, Diary-EEG Average), and save as .nii
%9. Perform correlations with summed connectivity to areas of a priori EEG-fMRI activation 

%Dependencies:

%nii_tool.m - a function for reading/writing .nii files (should be saved in /data/SOFTWARE)
%downloaded from: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/42997/versions/80/previews/nii_tool.m/index.html?access_key=

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

topdir = '/data/LGS-PP/MATLAB';
subject_list = load('/data/LGS_SUBJECTS.txt');
conndir = '/data/LGS-PP/VTA_CONNECTOMES';

addpath('/data/SOFTWARE'); %load in dependencies

alllgs_14subj_19vta_wdata_left_norm = zeros(19,429,14); %empty matrix to save mu-weighted data for LEFT vtas. Dimensions = 19 ESTEL patients x 429 regions x 14 disease-matched LGS subjects
alllgs_14subj_19vta_wdata_right_norm = zeros(19,429,14); %empty matrix to save mu-weighted data for RIGHT vtas. Dimensions = 19 ESTEL patients x 429 regions x 14 disease-matched LGS subjects

for s = 1:length(subject_list); %for each subject 
    
    for v = 1:19; %for each vta
        
         %load in connectome matrix for subject s, vta v
         subj=subject_list(s);
         
         datadir=fullfile(conndir, num2str(subj));
         data=load(fullfile(datadir, append('connectome_sift2_vta_',num2str(v),'.csv')));
         
         %remove first row/column (containing 'unassigned' streamlines)
         data = data(2:end, 2:end);
         
         %keep only columns 430:431 (corresponding to left vta and right vta, respectively)
         data = data(:,430:431);
         
         %remove vtas-as-targets connectivity values (430:431)
         data = data(1:429,:);
         
         %multiply each matrix element by coefficient 'mu' (output from sift2)
         mu = load(fullfile(datadir, 'sift2_mu.txt'));
         wdata = mu.*data;
         
         %norm by total connectivity strength for each vta
         wdata_norm = wdata./sum(wdata);
                  
         alllgs_14subj_19vta_wdata_left_norm(v,:,s) = wdata_norm(:,1)';
         alllgs_14subj_19vta_wdata_right_norm(v,:,s) = wdata_norm(:,2)'; 
         
    end
    
end

%average across 14 disease-matched lgs subjects (retaining per-vta values)
alllgs_subjavg_19vta_wdata_left_norm=mean(alllgs_14subj_19vta_wdata_left_norm,3); %3rd dimension = 14 disease-matched lgs subjects
alllgs_subjavg_19vta_wdata_right_norm=mean(alllgs_14subj_19vta_wdata_right_norm,3); %3rd dimension = 14 disease-matched lgs subjects

%sum together left and right vta lgs-subject-averaged matrices
alllgs_subjavg_19vta_wdata_bilatsum_norm=alllgs_subjavg_19vta_wdata_left_norm + alllgs_subjavg_19vta_wdata_right_norm;

%average across 19 vtas (yielding a single 429-length vector - one value per region)
alllgs_subjavg_vtaavg_wdata_bilatsum_norm=mean(alllgs_subjavg_19vta_wdata_bilatsum_norm);

%for display purposes, square-root normalise the lgs-subject-averaged & vta-averaged vector - to "compress" the small number of high values and "spread out" the
%larger number of low values
alllgs_subjavg_vtaavg_wdata_bilatsum_norm_sqrt=sqrt(alllgs_subjavg_vtaavg_wdata_bilatsum_norm);

%fill brain atlas with lgs-subject-averaged, vta-averaged, and sqrt-normalised values 
atlas = nii_tool('load', fullfile(topdir, 'MNI6THGEN_1MM_HCPMMP1_TIAN_CERES_AAN_1-429.nii')); %load in atlas
alllgs_subjavg_vtaavg_vol=zeros(size(atlas.img)); %blank matrix equal size to atlas, for filling with results

%fill with results
for iroi= 1 : 429;    
        alllgs_subjavg_vtaavg_vol(find(atlas.img==iroi))=alllgs_subjavg_vtaavg_wdata_bilatsum_norm_sqrt(iroi);
end

%save as nii
atlas_filled1=atlas;
atlas_filled1.img=alllgs_subjavg_vtaavg_vol;
nii_tool('save', atlas_filled1, 'alllgs_subjavg_vtaavg_wdata_bilatsum_norm_sqrt.nii');

%okay, now back to per-bilateral-vta, lgs-subject-averaged data. first exclude vtas for patients with missing eeg data (ESTEL patients 9 & 18)
alllgs_subjavg_19vta_wdata_bilatsum_norm_exc=alllgs_subjavg_19vta_wdata_bilatsum_norm;
alllgs_subjavg_19vta_wdata_bilatsum_norm_exc([9,18],:)=[];

%find mean and sqrt of above data (for later use re: connectivity to areas of a priori eeg-fmri activation)
alllgs_subjavg_vtaag_wdata_bilatsum_norm_exc=mean(alllgs_subjavg_19vta_wdata_bilatsum_norm_exc);
alllgs_subjavg_vtaag_wdata_bilatsum_norm_exc_sqrt=sqrt(alllgs_subjavg_vtaag_wdata_bilatsum_norm_exc);

%load in ESTEL patient clinical outcomes = each formatted as a column of percentages (% seizure reduction relative to baseline)
load(fullfile(topdir, 'lgs19_relreduc_diary_end3mths.mat')); %Diary data - 19 subjects - variable named "diary"
load(fullfile(topdir, 'lgs17_relreduc_EEG_end3mths.mat')); %EEG - 17 subjects (no 9 & 18) - variable named "EEG"
load(fullfile(topdir, 'lgs17_relreduc_EEGDiaryAvg_end3mths.mat')); %Diary-EEG Average - 17 subjects (no 9 & 18) - variable named "EEGDiaryAvg"

%correlate per-subject clinical outcomes with per-bilateral-vta lgs-subject-averaged connectivity strength
%DIARY OUTCOME !!!!
diary_conn_clin_corr=zeros(1,429);
for iroi=1:429;   
    conn=alllgs_subjavg_19vta_wdata_bilatsum_norm(:,iroi);
    RHO = corr(conn, diary, 'type', 'Spearman');
    diary_conn_clin_corr(iroi) = RHO;
end

%save rho values as .nii
rmap_szreduc_diary_vol=zeros(size(atlas.img)); %blank matrix equal size to atlas, for filling with results
for iroi=1:429;    
        rmap_szreduc_diary_vol(find(atlas.img==iroi))=diary_conn_clin_corr(iroi);
end

atlas_filled2=atlas;
atlas_filled2.img=rmap_szreduc_diary_vol;
nii_tool('save', atlas_filled2, 'rmap_spearman_diaryend3mths_wdata_bilatsum_norm.nii');

%correlate per-subj clinical outcomes with per-bilateral-vta lgs-subject-averaged connectivity strength
%EEG OUTCOME !!!!
eeg_conn_clin_corr=zeros(1,429);
for iroi=1:429;    
    conn=alllgs_subjavg_19vta_wdata_bilatsum_norm_exc(:,iroi); %note use of the *exc vector - patients with no eeg data excluded 
    RHO = corr(conn, EEG, 'type', 'Spearman');
    eeg_conn_clin_corr(iroi) = RHO;    
end

%save rho values as .nii
rmap_szreduc_eeg_vol=zeros(size(atlas.img)); %blank matrix equal size to atlas, for filling with results
for iroi=1:429;    
        rmap_szreduc_eeg_vol(find(atlas.img==iroi))=eeg_conn_clin_corr(iroi);
end

atlas_filled3=atlas;
atlas_filled3.img=rmap_szreduc_eeg_vol;
nii_tool('save', atlas_filled3, 'rmap_spearman_eegend3mths_17subj_wdata_bilatsum_norm.nii');

%correlate per-subj clinical outcomes with per-bilateral-vta lgs-subject-averaged connectivity strength
%DIARY-EEG AVG OUTCOME !!!!
eegdiaryavg_conn_clin_corr=zeros(1,429);
for iroi=1:429;   
    conn=alllgs_subjavg_19vta_wdata_bilatsum_norm_exc(:,iroi); %note use of the *exc vector - patients with no eeg data excluded 
    RHO = corr(conn, EEGDiaryAvg, 'type', 'Spearman');
    eegdiaryavg_conn_clin_corr(iroi) = RHO;    
end

%save rho values as .nii
rmap_szreduc_eegdiaryavg_vol=zeros(size(atlas.img)); %blank matrix equal size to atlas, for filling with results
for iroi=1:429;    
        rmap_szreduc_eegdiaryavg_vol(find(atlas.img==iroi))=eegdiaryavg_conn_clin_corr(iroi);
end

atlas_filled4=atlas;
atlas_filled4.img=rmap_szreduc_eegdiaryavg_vol;
nii_tool('save', atlas_filled4, 'rmap_spearman_eegdiaryavgend3mths_17subj_wdata_bilatsum_norm.nii');

%%%CORRELATION WITH AREAS OF A PRIORI EEG-fMRI GPFA-RELATED ACTIVATION

%load original and GPFA-masked versions of atlas
origatlas = nii_tool('load', fullfile(topdir, 'MNI6THGEN_1MM_HCPMMP1_TIAN_CERES_AAN_1-429.nii')); %load in atlas
pfaatlas = nii_tool('load', fullfile(topdir, 'MNI6THGEN_1MM_HCPMMP1_TIAN_CERES_AAN_1-429_GPFAMASKED.nii')); %load in atlas which has been masked with a binary verison of group-level eeg-fmri activation map 

%only keep atlas parcels where >50% of voxels remain after
%masking with the a priori EEG-fMRI activation map of generalised
%paroxysmal fast activity (GPFA) - as displayed in Figure 1C in paper

%to achieve this, find how many voxels contained within the original and
%gpfa-masked atlases, for each parcel. then set parcels to 0 if they contain fewer
%than 50% of original voxel number (before gpfa masking). 
origvoxno = [];
pfavoxno = [];
pfakeepatlas = [];
for iroi = 1:429;
    origvoxno(iroi) = length(find(origatlas.img == iroi));
    pfavoxno(iroi) = length(find(pfaatlas.img == iroi));
    if pfavoxno(iroi) >= origvoxno(iroi)*0.5; pfakeepatlas(iroi) = 1; else pfakeepatlas(iroi) = 0; end
end

%okay, now perform correlations with Diary-EEG average outcome: find parcels with rho>0 (0 = not >0, 1 = >0)
rpos_eegdiaryavg=eegdiaryavg_conn_clin_corr(1,:)>0;

%multiply by pfa keep index calculated above
rpos_pfakeep_eegdiaryavg=rpos_eegdiaryavg.*pfakeepatlas;

%get columns of connectivity matrix corresponding to rpos_pfakeep_diary
alllgs_subjavg_19vta_exc_eegdiaryavg_pfakeep=alllgs_subjavg_19vta_wdata_bilatsum_norm_exc(:,find(rpos_pfakeep_eegdiaryavg==1));

%sqrt-ed subject-average connectivity at each of the "kept" columns (for display purposes)
eediaryavg_rpos_pfakeep_connsubjavg_sqrt=rpos_pfakeep_eegdiaryavg.*alllgs_subjavg_vtaag_wdata_bilatsum_norm_exc_sqrt;

%save gpfa-retained connectivity values as .nii
eediaryavg_rpos_pfakeep_connsubjavg_sqrt_vol=zeros(size(atlas.img)); %blank matrix equal size to atlas, for filling with results
for iroi=1:429;    
        eediaryavg_rpos_pfakeep_connsubjavg_sqrt_vol(find(atlas.img==iroi))=eediaryavg_rpos_pfakeep_connsubjavg_sqrt(iroi);
end

atlas_filled5=atlas;
atlas_filled5.img=eediaryavg_rpos_pfakeep_connsubjavg_sqrt_vol;
nii_tool('save', atlas_filled5, 'alllgs_subjavg_vtaavg_wdata_bilatsum_norm_sqrt_17subj_eegdiaryavgend3mths_GPFAMASKED.nii');

%sum connectivity across columns (i.e., for each ESTEL patient, total
%connectivity to the gpfa-retained parcels)
eegdiaryavg_connsum_rpos_pfakeep=sum(alllgs_subjavg_19vta_exc_eegdiaryavg_pfakeep,2);

%lastly, perform one-tailed Spearman correlation between EEGDiaryAvg outcome and total connectivity to
%gpfa-retained areas
[r,p] = corr(eegdiaryavg_connsum_rpos_pfakeep, EEGDiaryAvg, 'type', 'Spearman', 'tail', 'right');