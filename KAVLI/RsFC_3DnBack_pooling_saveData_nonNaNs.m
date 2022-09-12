addpath(genpath('/space/syn02/1/data/MMILDB/mmildev/dhagler/mmil/MMPS/matlab/abcd'));
addpath(genpath('/home/abcdproj3/test/tfMRI_corr'));
addpath(['/home/abcdproj3/test/tfMRI_corr/MMPS_nonsvn']);
addpath(genpath('/home/wez025/matlab/ABCD_code/MultivariatePrediction'));

outputdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data/nonNaNs';
datadir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% final subjects that passed QC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/final_subjects.mat');
final_subjects = subjects_use;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(datadir, 'nback_data.mat'));
load(fullfile(datadir, 'nback_subjid.mat'));
nback_final_subject
tic
nonnans_idx = []
for idx = 1:size(nBack_data,1)
    tt = nBack_data{idx};
    nansums = sum(isnan(tt(:)));
    if nansums == 0 
        % sprintf('NaNs=%d, idx = %d\n', nansums, idx)
        nonnans_idx = [nonnans_idx, idx];
    end
end
toc

nback_final_subject = nback_final_subject(setdiff(1:length(nback_final_subject), nans_idx));
nBack_data = nBack_data(setdiff(1:length(nback_final_subject), nans_idx));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% rsfc - subject matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsfc = load('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/subjinfo.mat');
rsfc = rsfc.subjinfo;
rsfc_subjid_bl = rsfc.mri_info_visitid;

rsfmri = load('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/TSmat_pass.mat');
rsfmri = rsfmri.TSmat_pass;

idx = cell2mat(cellfun(@isempty, rsfmri, 'UniformOutput',false));


%{
 tic
save(fullfile(outputdir, 'rsfc_timeseries_data.mat'), 'rsfc_timeseries', '-v7.3');
toc

label = load('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/label.mat');
label = label.label;
label = label(ib);
save(fullfile(outputdir, 'label.mat'), 'label');

save(fullfile(outputdir, 'rsfc_subjid.mat'), 'rsfc_final_subjects');
 
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% nBack - convert 1D to 3D data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% [c ia ib] = intersect(final_subjects, nback.mri_info_visitid, 'stable');
% nback_final_subject = nback.mri_info_visitid(ib);
% nback = nback_dat(ib,:);

% fileID = fopen(fullfile('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data', 'test_idx.txt'),'r');
% formatSpec = '%d'
% test_idx = fscanf(fileID,formatSpec);
% size(test_idx)
% fclose(fileID);

% nBack_test = nback(test_idx,:);
% test_subj = final_subjects(test_idx);
% [c ia ib] = intersect('G031_INVM3E8EU28_baseline', test_subj, 'stable');
% test_subj = test_subj(595); %% 'G031_INVM3E8EU28_baseline'


nback_dat = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/tfmri/volmat_nBack_2_2_back_vs_0_back.mat');
nback_dat = nback_dat.volmat;

nback = load('/home/wez025/matlab/ABCD_code/KAVLI/data/vol_info_voxelwise_tfmri.mat');
nback = nback.subjinfo;
vol_mask_sub = nback.vol_mask_sub;
nback_subjid = nback.mri_info_visitid;

idx = find(ismember(nback.eventvec, 'baseline'));
nback_dat_bl = nback_dat(idx,:);
nback_subjid_bl = nback_subjid(idx);

ivec = find(isfinite(sum(nback_dat_bl, 2)));
nback_dat_bl = nback_dat_bl(ivec,:);

nback_subjid_bl = nback_subjid_bl(ivec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% dti - nonNaN baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dti = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/dmri/volinfo_dmri_VZ.mat');
dti = dti.dmri;
dti_subjid = dti.mri_info_visitid;
idx = find(ismember(dti.eventvec, 'baseline'));
dti_subjid_bl = dti_subjid(idx);

FA = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/dmri/volmat_FA.mat');
FA = FA.volmat;
FA_bl = FA(idx,:);
ivec_fa = find(isfinite(sum(FA_bl, 2)));
FA_bl = FA_bl(ivec_fa,:);

MD = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/dmri/volmat_MD.mat');
MD = MD.volmat;
MD_bl = MD(idx,:);
ivec_md = find(isfinite(sum(MD_bl, 2)));
MD_bl = MD_bl(ivec_fa,:);

dti_subjid_bl = dti_subjid_bl(ivec_md);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% dti - nonNaN baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smri = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/smri/volinfo_smri_VZ.mat');
smri = smri.smri;
smri_subjid = smri.mri_info_visitid;
for diri = 1:length(smri.mri_info_visitid)
    tmp = regexp(smri.mri_info_visitid{diri}, '^(?<site>\w+)_(?<SubjID>\w+)_(?<eventname>\w+)', 'names');
    smri.eventvec{diri} = tmp.eventname;
end
idx = find(ismember(smri.eventvec, 'baseline'));
smri_subjid_bl = smri_subjid(idx);

t1w = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/smri/volmat_nu.mat');
t1w = t1w.volmat;
t1w_bl = t1w(idx,:);
ivec = find(isfinite(sum(t1w_bl,2)));
t1w_bl = t1w_bl(ivec,:);
smri_subjid_bl = smri_subjid_bl(ivec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% getting overlapping subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

final_subjects = intersect(rsfc_subjid_bl, intersect(nback_subjid_bl, intersect(smri_subjid_bl, dti_subjid_bl,'stable')));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% nBack - convert 1D to 3D data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nBack_data = NaN([length(final_subjects), 100, 100, 130]);
nBack_data = cell([length(final_subjects), 1]);
tic
for idx = 1:length(final_subjects)
    subject_vector = nback_dat(ib(idx),:); 
    data_mask = subsample_volume(nback.vol_mask_aseg);
    nBack_data{idx} = fullvol(subject_vector, data_mask);
end
toc

tic
save(fullfile(outputdir, 'nback_data.mat'), 'nBack_data', '-v7.3');
toc

save(fullfile(outputdir, 'nback_data_subjid.mat'), 'nback_final_subject');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FA & MD - convert 1D to 3D data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dti = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/dmri/volinfo_dmri_VZ.mat');
dti = dti.dmri;
[c ia ib] = intersect(final_subjects, dti.mri_info_visitid, 'stable');
dti_final_subjects = dti.mri_info_visitid(ib);

FA = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/dmri/volmat_FA.mat');
FA = FA.volmat;
FA_data = cell([length(final_subjects) 1]);
tic
for idx = 1:length(final_subjects)
    subject_vector = FA(ib(idx),:); 
    data_mask = subsample_volume(dti.vol_mask_aseg);
    FA_data{idx} = fullvol(subject_vector, data_mask);
end
toc

tic
save(fullfile(outputdir, 'FA_data.mat'), 'FA_data', '-v7.3');
toc
save(fullfile(outputdir, 'FA_data_subjid.mat'), 'dti_final_subjects');


MD = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/dmri/volmat_MD.mat');
MD = MD.volmat;

MD_data = cell([length(final_subjects) 1]);
tic
for idx = 1:length(final_subjects)
    subject_vector = MD(ib(idx),:); 
    data_mask = subsample_volume(dti.vol_mask_aseg);
    MD_data{idx} = fullvol(subject_vector, data_mask);
end
toc

tic
save(fullfile(outputdir, 'MD_data.mat'), 'MD_data', '-v7.3');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% T1w, T2w - convert 1D to 3D data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smri = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/smri/volinfo_smri_VZ.mat');
smri = smri.smri;
[c ia ib] = intersect(final_subjects, smri.mri_info_visitid, 'stable');
smri_final_subjects = smri.mri_info_visitid(ib);

t1w = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/smri/volmat_nu.mat');
t1w = t1w.volmat;
t1w_data = cell([length(final_subjects) 1]);
tic
for idx = 1:length(final_subjects)
    subject_vector = t1w(ib(idx),:); 
    data_mask = subsample_volume(smri.vol_mask_aseg);
    t1w_data{idx} = fullvol(subject_vector, data_mask);
end
toc

tic
save(fullfile(outputdir, 't1w_data.mat'), 't1w_data', '-v7.3');
toc

save(fullfile(outputdir, 't1w_t2w_subjid.mat'), 'smri_final_subjects');


t2w = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/smri/volmat_T2.mat');
t2w = t2w.volmat;
t2w_data = cell([length(final_subjects) 1]);
tic
for idx = 1:length(final_subjects)
    subject_vector = t2w(ib(idx),:); 
    data_mask = subsample_volume(smri.vol_mask_aseg);
    t2w_data{idx} = fullvol(subject_vector, data_mask);
end
toc

tic
save(fullfile(outputdir, 't2w_data.mat'), 't2w_data', '-v7.3');
toc


JA = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/smri/volmat_JA.mat');
JA = JA.volmat;
JA_data = cell([length(final_subjects) 1]);
tic
for idx = 1:length(final_subjects)
    subject_vector = JA(ib(idx),:); 
    data_mask = subsample_volume(smri.vol_mask_aseg);
    JA_data{idx} = fullvol(subject_vector, data_mask);
end
toc

tic
save(fullfile(outputdir, 'JA_data.mat'), 'JA_data', '-v7.3');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% sanity check - checking all modalities have the same subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(outputdir, 't1w_t2w_subjid.mat'));
load(fullfile(outputdir, 'FA_MD_subjid.mat'));
load(fullfile(outputdir, 'nback_subjid.mat'));
load(fullfile(outputdir, 'rsfc_subjid.mat'));

[c ia ib] = intersect(rsfc_final_subjects, intersect(smri_final_subjects, intersect(dti_final_subjects, nback_final_subject, 'stable'), 'stable'), 'stable');


outputdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data';
load(fullfile(outputdir, 'nback_data.mat'));
t1 = nBack_data{1};
t2 = nBack_data{2};
t3 = nBack_data{3};

nback_subj = load(fullfile(outputdir, 'nback_subjid.mat'));
nback_subj = nback_subj.nback_final_subject;
[c ia ib] = intersect('S014_INVXHUFCVMW_baseline', nback_subj);