addpath(genpath('/space/syn02/1/data/MMILDB/mmildev/dhagler/mmil/MMPS/matlab/abcd'));
addpath(genpath('/home/abcdproj3/test/tfMRI_corr'));
addpath(['/home/abcdproj3/test/tfMRI_corr/MMPS_nonsvn']);
addpath(genpath('/home/wez025/matlab/ABCD_code/MultivariatePrediction'));

outputdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Baseline';

testdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data';
subjid = load(fullfile(testdir, 'nback_subjid.mat'));
subjid = subjid.nback_final_subject;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% get age and family ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = load('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/subjinfo.mat');
tmp = tmp.subjinfo;
[c ia ib] = intersect(subjid, tmp.mri_info_visitid, 'stable');
family_id = tmp.datamat_crossvalparam(ib);
age = tmp.datamat_covars(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% load nBack ROI level data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nback_ROI = sprintf('%s/nback_ROI.csv', '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/');
nback_ROI = readtable(nback_ROI, 'TreatAsEmpty','NA');

nback_ROI_subjid = nback_ROI.src_subject_id;
nback_ROI_data = nback_ROI(:,6:size(nback_ROI,2));
nback_ROI_data = cell2mat(table2cell(nback_ROI_data));

nback_subj = load('/home/wez025/matlab/ABCD_code/KAVLI/data/vol_info_voxelwise_tfmri.mat');
nback_subj = nback_subj.subjinfo;
baseline_idx = ismember(nback_subj.eventvec, 'baseline');
nback_subj_subjid = nback_subj.subjidvec(baseline_idx);
nback_subj_visitid = nback_subj.mri_info_visitid(baseline_idx);

[c ia ib] = intersect(nback_ROI_subjid, nback_subj_subjid, 'stable');
nback_ROI_mri_info_visitid = nback_subj_visitid(ib);
nback_ROI_subjid = nback_ROI_subjid(ia);
nback_ROI_data = nback_ROI_data(ia,:);

%%%%%%%
testdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/Data';
final_subjects = load(fullfile(testdir, 'nback_subjid.mat'));
final_subjects = final_subjects.nback_final_subject;

fileID = fopen(fullfile(testdir, 'test_idx.txt'),'r');
formatSpec = '%d'
test_idx = fscanf(fileID,formatSpec);
size(test_idx)
fclose(fileID);

train_idx = setdiff(1:length(subjid), test_idx);

test_subj = subjid(test_idx);
train_subj = subjid(train_idx);

label = load(fullfile(testdir, 'label.mat'));
label = label.label;
test_label = label(test_idx);
train_label = label(train_idx);
test_age = age(test_idx);
train_age = age(train_idx);
test_family_id = family_id(test_idx);
train_family_id = family_id(train_idx);

[c ia ib] = intersect(test_subj, nback_ROI_mri_info_visitid, 'stable');
test_subj = test_subj(ia);
test_label = test_label(ia);
test_age = test_age(ia,:);
test_family_id = test_family_id(ia);

test_data = nback_ROI_data(ib,:);

[c ia ib] = intersect(train_subj, nback_ROI_mri_info_visitid, 'stable');
train_subj = train_subj(ia);
train_label = train_label(ia);
train_age = train_age(ia,:);
train_family_id = train_family_id(ia);

train_data = nback_ROI_data(ib,:);


defvec = isfinite(sum(train_label,2)+sum(train_data,2)); ivec = find(defvec); %% finding subjects with complete data for all matrices
train_label = train_label(ivec);
train_data = train_data(ivec,:);
train_family_id = train_family_id(ivec);
train_subj = train_subj(ivec);
train_age = train_age(ivec);

training_output = ABCD_prediction_svd_covars_5foldcrossvalwithin(train_label, train_data, train_family_id, train_subj, 'agevar', train_age);
kval = mean(training_output.maxK_mat);

defvec = isfinite(sum(test_label,2)+sum(test_data,2)); ivec = find(defvec); %% finding subjects with complete data for all matrices
test_label = test_label(ivec);
test_data = test_data(ivec,:);
test_family_id = test_family_id(ivec);
test_subj = test_subj(ivec);
test_age = test_age(ivec);

[beta_hat beta_se zmat logpmat sig2tvec sig2mat binvec] = SSE_fit(valvec_disc_subfold,iid_subfold,fam_id_subfold,agevec_subfold,valmat_disc_subfold);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% get age and family ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = load('/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/subjinfo.mat');
tmp = tmp.subjinfo;
[c ia ib] = intersect(subjid, tmp.mri_info_visitid, 'stable');
family_id = tmp.datamat_crossvalparam(ib);
age = tmp.datamat_covars(:,1);



[c ia ib] = intersect(subjid, nback_ROI_mri_info_visitid, 'stable');

subjid = subjid(ia);
family_id = family_id(ia);
age = age(ia);
nback_ROI_data = nback_ROI_data(ib,:);
nback_ROI_mri_info_visitid = nback_ROI_mri_info_visitid(ib);
nback_ROI_subjid = nback_ROI_subjid(ib);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% load vertexwise data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertexdir = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/vertexwise/tfmri';
imgname = 'nBack_2_back_vs_0_back';

load('/home/wez025/matlab/ABCD_release_3.0/SurfeView_surfs.mat');
icnum = 4; icnvert = size(icsurfs{icnum}.vertices,1); % Icosahedral order to use for stats
tmp = regexp(char(imgname), 'w*_', 'split', 'once'); modname = tmp{1}; condname = tmp{2};
hemistrings = {'lh','rh'};

betamat = [];
for hemii = 1:2
    hemi = hemistrings{hemii};
    fname = sprintf('%s/%s_1_%s-%s.mat',vertexdir, modname, condname,hemi); % Should save these pre-truncated to specified icnum
    tmp = load(fname);
    betamat = cat(2,betamat,tmp.betamat(:,1:size(icsurfs{icnum}.vertices,1)));
end

vertex_sub = load(fullfile(vertexdir, 'vol_info.mat'));
vertex_sub = vertex_sub.visitidvec;

[c ia ib] = intersect(subjid, vertex_sub, 'stable');
betamat = betamat(ib,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% read testset id
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen(fullfile(testdir, 'test_idx.txt'),'r');
formatSpec = '%d'
test_idx = fscanf(fileID,formatSpec);
size(test_idx)
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% matching training and testing set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
train_idx = setdiff(1:length(subjid), test_idx);

test_subj = subjid(test_idx);
train_subj = subjid(train_idx);

label = load(fullfile(testdir, 'label.mat'));
label = label.label;

test_label = label(test_idx);
train_label = label(train_idx);

test_nback = betamat(test_idx,:);
train_nback = betamat(train_idx,:);

test_age = age(test_idx);
train_age = age(train_idx);

test_family_id = family_id(test_idx);
train_family_id = family_id(train_idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% doing prediction!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defvec = isfinite(sum(train_label,2)+sum(train_nback,2)); ivec = find(defvec); %% finding subjects with complete data for all matrices
train_label = train_label(ivec);
train_nback = train_nback(ivec,:);
train_family_id = train_family_id(ivec);
train_subj = train_subj(ivec);
train_age = train_age(ivec);

output = ABCD_prediction_svd_covars_5foldcrossvalwithin(train_label, train_nback, train_family_id, train_subj, 'agevar', train_age);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% load voxelwise data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dat = load(fullfile(testdir, 'nback_data.mat'));
nback_dat = load('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/imaging_concat/voxelwise/ABCD2_cor10/tfmri/volmat_nBack_1_2_back_vs_0_back.mat');
nback_dat = nback_dat.volmat;

nback_subj = load('/home/wez025/matlab/ABCD_code/KAVLI/data/vol_info_voxelwise_tfmri.mat');
nback_subj = nback_subj.subjinfo;
% vol_mask_sub = nback_subj.vol_mask_sub;
nback_subjid = nback_subj.mri_info_visitid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% read testset id
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(fullfile(testdir, 'test_idx.txt'),'r');
formatSpec = '%d'
test_idx = fscanf(fileID,formatSpec);
size(test_idx)
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% subject matching  & train and test set separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

train_idx = setdiff(1:length(subjid), test_idx);

test_subj = subjid(test_idx);
train_subj = subjid(train_idx);

label = load(fullfile(testdir, 'label.mat'));
label = label.label;

test_label = label(test_idx);
train_label = label(train_idx);

[c ia ib] = intersect(test_subj, nback_subjid, 'stable');
test_nback = nback_dat(ib,:);

[c ia ib] = intersect(train_subj, nback_subjid, 'stable');
train_nback = nback_dat(ib,:);

test_age = age(test_idx);
train_age = age(train_idx);

test_family_id = family_id(test_idx);
train_family_id = family_id(train_idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% doing prediction!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defvec = isfinite(sum(train_label,2)+sum(train_nback,2)); ivec = find(defvec); %% finding subjects with complete data for all matrices
train_label = train_label(ivec);
train_nback = train_nback(ivec,:);
train_family_id = train_family_id(ivec);
train_subj = train_subj(ivec);
train_age = train_age(ivec);
% 
% output = ABCD_prediction_svd_covars_5foldcrossvalwithin(train_label, train_nback, train_family_id, train_subj, 'agevar', train_age);


rng(1)
Lambda = [0.0001:0.05:1];

train_nback = train_nback'; 
% CVMdl = fitrlinear(train_nback, train_label,'ObservationsIn','columns','KFold',5,'Lambda',Lambda,...
%     'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
% mse = kfoldLoss(CVMdl);
YHat = kfoldPredict(CVMdl);

[Mdl,FitInfo] = fitrlinear(train_nback,train_label);

hyperopts = struct('AcquisitionFunctionName','expected-improvement-plus');
[Mdl,FitInfo,HyperparameterOptimizationResults] = fitrlinear(train_nback,train_label,'ObservationsIn', 'columns', ...
    'OptimizeHyperparameters','all',...
    'HyperparameterOptimizationOptions',hyperopts)
YHat = predict(Mdl, test_nback);
corr(test_label, YHat)

%%% random forest with PCA dimension reduction
rng(1)
[coeff,score,latent,tsquared,explained,mu] = pca(train_nback, 'Economy', true, 'NumComponents',1);
[Mdl,FitInfo] = fitrlinear(score,train_label);

corr(score*Mdl.Beta, train_label)

scoreTest95 = (test_nback-mu)*coeff(:,1);
YTest_predicted = predict(Mdl,scoreTest95);
mdl = fitlm(score,train_label);

RF_Mdl = TreeBagger(200,train_nback,train_label, 'Method','regression');



lambda_final = Lambda(find(min(mse)));
lasso_mod = fitrlinear(train_nback, train_label,'ObservationsIn','columns','Lambda',lambda_final,...
'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
final_betas = lasso_mod.Beta;

%%%%% test set
defvec = isfinite(sum(test_label,2)+sum(test_nback,2)+sum(test_age,2)); ivec = find(defvec); %% finding subjects with complete data for all matrices
test_label = test_label(ivec);
test_nback = test_nback(ivec,:);

test_label_pred = test_nback*final_betas;
corr(test_label_pred, test_label)
