addpath(genpath('/space/syn02/1/data/MMILDB/mmildev/dhagler/mmil/MMPS/matlab/abcd'));
addpath(genpath('/home/abcdproj3/test/tfMRI_corr'));
addpath(['/home/abcdproj3/test/tfMRI_corr/MMPS_nonsvn']);
addpath(genpath('/home/wez025/matlab/ABCD_code/MultivariatePrediction'));

outputdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% final subject list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/final_subjects.mat');
final_subjects = subjects_use;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% rsfc - subject matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsfc = load('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/subjinfo.mat');
rsfc = rsfc.subjinfo;
[c ia ib] = intersect(final_subjects, rsfc.mri_info_visitid, 'stable');
rsfc_final_subjects = rsfc.mri_info_visitid(ib);

rsfmri = load('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/TSmat_pass.mat');
rsfmri = rsfmri.TSmat_pass;
rsfc_timeseries = rsfmri(ib);
tic
save(fullfile(outputdir, 'rsfc_timeseries_data.mat'), 'rsfc_timeseries', '-v7.3');
toc

label = load('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/label.mat');
label = label.label;
label = label(ib);
save(fullfile(outputdir, 'label.mat'), 'label');

save(fullfile(outputdir, 'rsfc_subjid.mat'), 'rsfc_final_subjects');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% nBack - convert 1D to 3D data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nback = load('/home/wez025/matlab/ABCD_code/KAVLI/data/vol_info_voxelwise_tfmri.mat');
nback = nback.subjinfo;
vol_mask_sub = nback.vol_mask_sub;
nback_subjid = nback.mri_info_visitid;

[c ia ib] = intersect(final_subjects, nback.mri_info_visitid, 'stable');
nback_final_subject = nback.mri_info_visitid(ib);

nback_dat = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/tfmri/volmat_nBack_2_2_back_vs_0_back.mat');
nback_dat = nback_dat.volmat;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% cbcl label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cbcl = sprintf('%s/cbcl_baseline_4.0.csv', '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data/');
cbcl = readtable(cbcl, 'TreatAsEmpty','NA');

dti = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/dmri/volinfo_dmri_VZ.mat');
dti = dti.dmri;

bl_idx = find(ismember(dti.eventvec, 'baseline'));
dti_bl_mri_info_visitid = dti.mri_info_visitid(bl_idx);
dti_bl_subjid = dti.subjidvec(bl_idx);
dti_bl_subjid = strcat('NDAR_', dti_bl_subjid);

[c ia ib] = intersect(cbcl.src_subject_id, dti_bl_subjid, 'stable');
cbcl = cbcl(ia,:);
cbcl.mri_info_visitid = dti_bl_mri_info_visitid(ib)';


cbcl_totprob = cbcl.cbcl_scr_syn_totprob_r;
cbcl_externalprob = cbcl.cbcl_scr_syn_external_r;
cbcl_internalprob = cbcl.cbcl_scr_syn_internal_r;
cbcl_mri_info_visitid = cbcl.mri_info_visitid;

MH = struct;
MH.cbcl_totprob = cbcl_totprob;
MH.cbcl_externalprob = cbcl_externalprob;
MH.cbcl_internalprob = cbcl_internalprob;
MH.cbcl_mri_info_visitid = cbcl_mri_info_visitid;


save(fullfile('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data/', 'MH.mat'), 'MH');

load(fullfile('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data/', 'MH.mat'));
cbcl_totprob = MH.cbcl_totprob;
cbcl_subjid = MH.cbcl_mri_info_visitid;

subjid = load(fullfile('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data/', 'rsfc_subjid.mat'));
subjid = subjid.rsfc_final_subjects;

[c ia ib] = intersect(cbcl_subjid, subjid, 'stable');
cbcl_totprob = cbcl_totprob(ia);
cbcl_subjid = cbcl_subjid(ia);

[c ia ib] = intersect(cbcl_subjid, subjid, 'stable');
cbcl_totprob_label = NaN([size(subjid)]);
cbcl_totprob_label(ib) = cbcl_totprob;

save(fullfile('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data/', 'label_cbcl_totprob.mat'), 'cbcl_totprob_label');





