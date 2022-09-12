addpath(genpath('/space/syn02/1/data/MMILDB/mmildev/dhagler/mmil/MMPS/matlab/abcd'));
addpath(genpath('/home/abcdproj3/test/tfMRI_corr'));
addpath(['/home/abcdproj3/test/tfMRI_corr/MMPS_nonsvn']);
addpath(genpath('/home/wez025/matlab/ABCD_code/MultivariatePrediction'));

% For R4.0
% datadirs = {'/space/syn14/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn16/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn18/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn20/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn24/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn27/1/data/ABCD/DAL_ABCD_R4/proc_bold' ...
%             '/space/syn28/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn29/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn31/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn34/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn35/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn36/1/data/ABCD/DAL_ABCD_R4/proc_bold' ...
%              '/space/syn39/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn44/1/data/ABCD/DAL_ABCD_R4/proc_bold' '/space/syn52/1/data/ABCD/DAL_ABCD_R4/proc_bold' ...
%            };
% %datadirs = cellfun(@(x) sprintf('/space/%s/1/data/ABCD/DAL_ABCD%s/proc_bold',sprintf('syn%02d',x),sprintf('_syn%02d',x)),num2cell([1:50]'),'UniformOutput',false); % Look in every syn* file system
%
load(fullfile('/space/gwas-syn2/1/data/GWAS/ABCD/Images/NDA3.0/rsFC_R4/R4_newpreproc/', sprintf('Rs_roi_newpreproc_roinames.mat')));
load(fullfile('/space/gwas-syn2/1/data/GWAS/ABCD/Images/NDA3.0/rsFC_R4/R4_newpreproc/', sprintf('Release4.0_Rs_roi_subinfo_baseline.mat')));
SubjIDvec_inventory = subjinfo.SubjIDvec_inventory;
sitevec_inventory = subjinfo.sitevec_inventory;
datevec_inventory = subjinfo.datevec_inventory;
eventname_inventory = subjinfo.eventname_inventory;
VisitIDvec_inventory = subjinfo.VisitIDvec_inventory;
dirvec = subjinfo.dirvec_inventory;
ndirs = length(dirvec);

% Rs_run = cell([ndirs 1]); 
% nts_run = cell(size(Rs_run)); 
% meanmotion_run = cell(size(Rs_run));

% diristop = 1;
% % infix = 'corr_resBOLD';
% for diri = diristop:ndirs
%     if mod(diri,100) == 0
%         fprintf(1,'\n\n%d/%d (%s)\n\n',diri,ndirs,datestr(now));
%     end
%     indir = dirvec{diri};
%     try
%         data = abcd_load_rsfMRI_roi_data_VZ(indir);
%     catch
%         fprintf(1,'%s: file missing or corrupted\n',indir);
%     end
%     if ~isempty(data)
%         data_prep = abcd_prep_rsfMRI_roi_data_db(data);

%         %%% filtering ROIs: Gordon parcel and subcortical ROIs.
%         jvec_all = [1:length(roinames)];
%         jvec_gp = find(find_cell(strfind(roinames,'-gp'))); % Find Gordon parcels
%         jvec_dk = find(find_cell(strfind(roinames,'ctx-'))&~find_cell(strfind(roinames,'-gp'))); % Find FS Desikan parcels
%         jvec_subcort = [555:559,562:564,567,568,574:580,582:583];

%         jvec_use = roinames([jvec_gp, jvec_subcort]); % Plot GP parcels only

%         D = []; ntsum = []; meanmotion = [];
%         nscans = length(data_prep);
%         for r=1:nscans
%             ntpoints = length(data_prep(r).ind_valid);
%             if ntpoints>=20
%                 motion_radius = 50;
%                 motion_absflag = true;
%                 censor_thresh = 0.2;
%                 motion_nodyflag = false;
%                 [motion_stats motion_fd] = mmil_motion_stats(data_prep(r).motion_data, motion_radius, motion_absflag, censor_thresh, motion_nodyflag);
%                 motion = motion_stats.mean_motion;
%                 roi_use = find(ismember(data_prep(r).roinames, jvec_use));

%                 D = cat(1, D, data_prep(r).roi_data_censored(:,roi_use));
%                 ntsum = cat(1, ntsum, ntpoints);
%                 meanmotion = cat(1, meanmotion, mean(motion_fd));
%             else 
%                 continue;
%             end;
%         end;  
%         Rs_run{diri} = D;
%         nts_run{diri} = ntsum;
%         meanmotion_run{diri} = meanmotion;  
%     else
%             fprintf(1,'No rsFC data for %s\n', indir)
%     end
%     if mod(diri, 500) == 0
%         save(fullfile(outputdir, sprintf('Release4.0_newpreproc_Rs_roi_baseline.mat')), 'Rs_run', '-v7.3')
%         save(fullfile(outputdir, sprintf('Release4.0_newpreproc_Rs_roi_nts_baseline.mat')), 'nts_run')
%         save(fullfile(outputdir, sprintf('Release4.0_newpreproc_Rs_roi_meanmotion_baseline.mat')), 'meanmotion_run')
%     end
% end
% save(fullfile(outputdir, sprintf('Release4.0_newpreproc_Rs_roi_baseline.mat')), 'Rs_run', '-v7.3')
% save(fullfile(outputdir, sprintf('Release4.0_newpreproc_Rs_roi_nts_baseline.mat')), 'nts_run')
% save(fullfile(outputdir, sprintf('Release4.0_newpreproc_Rs_roi_meanmotion_baseline.mat')), 'meanmotion_run')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

datadir = '/space/gwas-syn2/1/data/GWAS/ABCD/Images/NDA3.0/rsFC_R4/R4_newpreproc/nback_rsfc/separateRuns/ntsbyroi';
tt = load(fullfile(datadir, 'Release4.0_newpreproc_Rs_roi_baseline.mat'));
tmp = tt.Rs_run;
meanmotion = load(sprintf('%s/Release4.0_newpreproc_Rs_roi_meanmotion_baseline.mat', datadir)); meanmotion = meanmotion.meanmotion_run;
nts = load(sprintf('%s/Release4.0_newpreproc_Rs_roi_nts_baseline.mat', datadir)); nts = nts.nts_run;

nts = cellfun(@sum, nts);
meanmotion = cellfun(@mean, meanmotion);

subjinfo = load(sprintf('%s/Release4.0_Rs_roi_subinfo_baseline.mat', '/space/gwas-syn2/1/data/GWAS/ABCD/Images/NDA3.0/rsFC_R4/R4_newpreproc/')); subjinfo = subjinfo.subjinfo;
VisitIDvec_inventory = subjinfo.VisitIDvec_inventory;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% concatenate the correlation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nrows, ncols] = cellfun(@size, tmp);

Rmat = NaN([size(tmp,1),352, 352]);
Zmat = NaN(size(Rmat));

for subj = 1:size(tmp,1)
    D = tmp{subj};
    if ~isempty(D)
        if size(D,1)<375
            continue;
        else
            R = corr(D);
            Z = atanh(R);
            Z(~isfinite(Z)) = 0;
            Rmat(subj,:,:) = R;
            Zmat(subj,:,:) = Z;
        end
    else
        fprintf('WARNING: no valid data\n');
    end
end

nonempty_idx = find(mean(mean(~isnan(Rmat),2),3)==1);
VisitIDvec_notnan = VisitIDvec_inventory(nonempty_idx);
Rmat_pass = Rmat(nonempty_idx,:,:);
Zmat_pass = Zmat(nonempty_idx,:,:);
nts = nts(nonempty_idx);
meanmotion = meanmotion(nonempty_idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% match and save the ROIxtimeseries data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
include_idx = [];
for i = 1:length(tmp)
    tt = tmp{i};
    if ~isempty(tt) & size(tt,1)>=375
        include_idx = [include_idx, i];
    end
end

TSmat_pass = tmp(include_idx);
nts = nts(include_idx);
meanmotion = meanmotion(include_idx);
VisitIDvec_notnan = VisitIDvec_inventory(include_idx);
sprintf('%d %d %d %d', size(TSmat_pass,1), size(nts,1), size(meanmotion,1), size(VisitIDvec_notnan,1))

% nan_idx = find(cell2mat(cellfun(@isempty, tmp, 'UniformOutput', false)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% behavioral phenotypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tabulateddir = '/home/wez025/matlab/ABCD_release_3.0/Released/Data'; %% behavioral measures and covariates

behavdat = sprintf('%s/behVars_release4.0.csv', tabulateddir);
behavior_mat = readtable(behavdat, 'TreatAsEmpty','NA');
colname_BEHAV = behavior_mat.Properties.VariableNames;
%%% create iids variable for matching
behavior_mat.iids = strcat(behavior_mat.src_subject_id,'_',behavior_mat.eventname);
%% sort and select baseline data only
baseline_ind = ismember(behavior_mat.eventname,'baseline_year_1_arm_1');
behavior_mat = behavior_mat(baseline_ind,:);

colilist = [37, 63]; %% 37=nihtbx_totalcomp_uncorrected, 63=ssrt
behavtasklist = colname_BEHAV(colilist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% qcflags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname_tabulated = sprintf('%s/qcflags_dof_nts_meanmotion_release4.0.csv', '/home/wez025/matlab/ABCD_release_3.0/Released/Data');
fprintf(1,'Reading tabulated imaging data from %s\n',fname_tabulated);
% opts = detectImportOptions(fname_tabulated, 'NumHeaderLines',1);
qctable_mat = readtable(fname_tabulated,'ReadVariableNames', true, 'TreatAsEmpty','NA');
ind = ismember(qctable_mat.eventname, 'baseline_year_1_arm_1');
ind = ismember(qctable_mat.mri_info_visitid, 'G010_INV4AY58X03_baseline');
qctable_mat.imgincl_nback_include(2747)
qctable_mat.imgincl_rsfmri_include(2747)


qctable_mat = qctable_mat(ind,:);
qctable_mat.iids_img = strcat(qctable_mat.src_subject_id, '_', qctable_mat.eventname);

load(fullfile('/space/gwas-syn2/1/data/GWAS/ABCD/Images/NDA3.0/taskrsfc_R4/', sprintf('Release4.0_taskresid_Rs_roi_subinfo_baseline.mat')));
VisitIDvec_inventory = subjinfo.VisitIDvec_inventory;
visitid_eventname = strcat(subjinfo.sitevec_inventory, '_', subjinfo.SubjIDvec_inventory, '_', subjinfo.eventname_inventory);

[c ia ib] = intersect(qctable_mat.mri_info_visitid, visitid_eventname,'stable');
qctable_mat = qctable_mat(ia,:);
qctable_mat.mri_info_visitid_eventname = qctable_mat.mri_info_visitid;
qctable_mat.mri_info_visitid = VisitIDvec_inventory(ib);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% covariates; contains INTERCEPT (imputed income)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fname_tabulated = sprintf('%s/ABCD_designmatrix.csv', tabulateddir); %% original income
fname_tabulated = sprintf('%s/ABCD_designmatrix_incomeimp_release4.0.csv', tabulateddir); %% original income
covars_mat = readtable(fname_tabulated,'ReadVariableNames', true);
ind = ismember(covars_mat.eventname, 'baseline_year_1_arm_1');
covars_mat = covars_mat(ind,:);
% iid_covars = covars.src_subject_id;
% eid_covars = covars.eventname;
covars_mat.iids_covars = strcat(covars_mat.src_subject_id, '_', covars_mat.eventname);
tt = covars_mat.Properties.VariableNames;

covar_ind = cell([4, 1]); 

covar_ind{1} = [16:51,59]; %% scanner and software variables
covar_ind{2} = [4,15:51,59]; %% scanner/software, age, sex
covar_ind{3} = [4:51,58,59]; %% scanner/software, age, sex, top 10 genetic PCs
covar_ind{4} = [4:59]; %% scanner/software, age, sex, top 10 genetic PCs, householdi income, parental education
covarnames = {'Scan', 'ageSexScan', 'ageSexScanPCs', 'SES'};
ncovars = size(covar_ind,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% matching subjects with complete behavior and covariates data, and have passed rsfc_qc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
behav = behavior_mat;
colname_BEHAV = behav.Properties.VariableNames; INVlist_meta = behav.iids;
datavec_BEHAV = table2cell(behav); 
% loading covariates
covars = covars_mat;
colname_COVARS = covars.Properties.VariableNames;
datamat_crossvalparam = covars.rel_family_id;  
INVlist_covars = covars.iids_covars;
covars = table2cell(covars); 
datamat_covars = NaN(size(covars,1), length(covar_ind{4})); 
for i = 1:length(covar_ind{4})
    ci = covar_ind{4}(i);
    datamat_covars(:,i) = cell2mat(covars(:,ci));
end
covars_name = colname_COVARS(covar_ind{4});

[dummy I_meta I_covars] = intersect(INVlist_meta,INVlist_covars,'stable');
INVlist_covars = INVlist_covars(I_covars);
datamat_covars = datamat_covars(I_covars,:);
datamat_crossvalparam = datamat_crossvalparam(I_covars);
INVlist_meta = INVlist_meta(I_meta);

ncols = length(colilist);
datamat_meta = NaN(length(INVlist_meta),ncols);
for ci = 1:ncols
    datamat_meta(:,ci) = cell2mat(datavec_BEHAV(I_meta,colilist(ci)));
end

%%%%%%%%%% loading qc data
qctable = qctable_mat;
colnames_qc = qctable.Properties.VariableNames; 
qc_flag = find(contains(colnames_qc, 'imgincl_rsfmri_include'));
qctable = table2cell(qctable); 
vertmatchcol = find(contains(colnames_qc,'mri_info_visitid')); %% subject IDs to match with imaging data, starts with "S30_..."
behmatchcol = find(contains(colnames_qc,'iids_img')); %% subject IDs to match with behavioral data, starts with "NDAR..."
qctable = qctable(find(cell2mat(qctable(:,qc_flag))== 1),:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% matching brain data with rsfc_qc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dummy I_vert, I_dof] = intersect(VisitIDvec_notnan, qctable(:,vertmatchcol), 'stable');
VisitIDvec_notnan = VisitIDvec_notnan(I_vert);
% Rmat_pass = Rmat_pass(I_vert,:,:);
% Zmat_pass = Zmat_pass(I_vert,:,:);
TSmat_pass = TSmat_pass(I_vert);
nts = nts(I_vert);
meanmotion = meanmotion(I_vert);
qctable = qctable(I_dof,:);

%%%%%%%%% matching brain and behavior
[dummy I_meta I_vert] = intersect(INVlist_meta, qctable(:,behmatchcol),'stable');

datamat_covars = datamat_covars(I_meta,:);
subjID = qctable(:,vertmatchcol);
subjID = subjID(I_vert);
datamat_crossvalparam = datamat_crossvalparam(I_meta);
datamat_meta = datamat_meta(I_meta,1);

% Rmat_pass = Rmat_pass(I_vert,:,:);
% Zmat_pass = Zmat_pass(I_vert,:,:);
TSmat_pass = TSmat_pass(I_vert);
nts = nts(I_vert);
meanmotion = meanmotion(I_vert);

subjinfo = struct();
subjinfo.subjID = subjID;
subjinfo.datamat_crossvalparam = datamat_crossvalparam;
subjinfo.datamat_covars = datamat_covars;
subjinfo.nts = nts;
subjinfo.meanmotion = meanmotion;

label = datamat_meta;

outputdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI';
if exist(outputdir)
    mkdir(outputdir)
end
save(fullfile(outputdir, 'label.mat'), 'label');
save(fullfile(outputdir, 'subjinfo.mat'), 'subjinfo');
save(fullfile(outputdir, 'TSmat_pass.mat'), 'TSmat_pass', '-v7.3');
% save(fullfile(outputdir, 'Rmat_pass.mat'), 'Rmat_pass', '-v7.3');
% save(fullfile(outputdir, 'Zmat_pass.mat'), 'Zmat_pass', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% get network IDs
gp = readtable('/home/cfan/codes/ABCD_analysis/gordon_parcels.csv');
jvec = [217:549,554:559,562:564,572:581]; % gp + subcort
lookup = table(jvec',roinames(jvec)', cat(1, gp.community, repmat({'Subcort'}, 19, 1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% task-state FC concatenation
%%%%%%%%%%%%%%%%%%%%% pass qc for that modality, save subjectID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taskrsfc_R4 = '/space/gwas-syn2/1/data/GWAS/ABCD/Images/NDA3.0/taskrsfc_R4/nback_rsfc/separate_runs/rawtimeseries/ntsbyroi'; %% taskresidualized FC

imgname = {'nBack', 'SST', 'MID'};

for taski = 2:3

    load(fullfile('/space/gwas-syn2/1/data/GWAS/ABCD/Images/NDA3.0/taskrsfc_R4/', sprintf('Release4.0_taskresid_Rs_roi_subinfo_baseline.mat')));
    VisitIDvec_inventory = strcat(subjinfo.sitevec_inventory, '_', subjinfo.SubjIDvec_inventory, '_', subjinfo.eventname_inventory);
    dirvec = subjinfo.dirvec_inventory;
    ndirs = length(dirvec);

    data = load(fullfile(taskrsfc_R4, sprintf('Release4.0_%s_rawtimeseries_taskresid_Rs_run_baseline_ntsbyroi.mat', imgname{taski})));
    tmp = data.Rs_run;
    nts = load(fullfile(taskrsfc_R4, sprintf('Release4.0_%s_rawtimeseries_taskresid_nts_baseline_ntsbyroi.mat', imgname{taski})));
    nts = nts.nts_run;
    meanmotion = load(fullfile(taskrsfc_R4, sprintf('Release4.0_%s_rawtimeseries_taskresid_meanmotion_baseline_ntsbyroi.mat', imgname{taski})));
    meanmotion = meanmotion.meanmotion_run;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% setting subjects not passing QC to NaN 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fname_tabulated = sprintf('%s/qcflags_dof_nts_meanmotion_release4.0.csv', '/home/wez025/matlab/ABCD_release_3.0/Released/Data');
    fprintf(1,'Reading tabulated imaging data from %s\n',fname_tabulated);
    % opts = detectImportOptions(fname_tabulated, 'NumHeaderLines',1);
    qctable_mat = readtable(fname_tabulated,'ReadVariableNames', true, 'TreatAsEmpty','NA');
    ind = ismember(qctable_mat.eventname, 'baseline_year_1_arm_1');
    qctable_mat = qctable_mat(ind,:);
    qctable_mat.iids_img = strcat(qctable_mat.src_subject_id, '_', qctable_mat.eventname);

    qctable = qctable_mat;
    colnames_qc = qctable.Properties.VariableNames; 
    if ismember(imgname{taski}, {'nBack'})
        qc_flag = find(contains(colnames_qc, 'imgincl_taskres_nback_include'));
    elseif ismember(imgname{taski}, {'SST'})
        qc_flag = find(contains(colnames_qc, 'imgincl_taskres_sst_include'));
    elseif ismember(imgname{taski}, {'MID'})
        qc_flag = find(contains(colnames_qc, 'imgincl_taskres_mid_include'));
    end
    qctable = table2cell(qctable); 
    vertmatchcol = find(contains(colnames_qc,'mri_info_visitid')); %% subject IDs to match with imaging data, starts with "S30_..."
    behmatchcol = find(contains(colnames_qc,'iids_img')); %% subject IDs to match with behavioral data, starts with "NDAR..."
    qctable = qctable(find(cell2mat(qctable(:,qc_flag))== 1),:); 
    
    [dummy I_vert, I_dof] = intersect(VisitIDvec_inventory, qctable(:,vertmatchcol), 'stable');
    VisitIDvec_inventory = VisitIDvec_inventory(I_vert);
    % Rmat_pass = Rmat_pass(I_vert,:,:);
    % Zmat_pass = Zmat_pass(I_vert,:,:);
    tmp = tmp(I_vert);
    nts = nts(I_vert);
    meanmotion = meanmotion(I_vert);
    qctable = qctable(I_dof,:);

    nonempty_idx = find(~cellfun(@isempty, tmp));
    VisitIDvec_inventory = VisitIDvec_inventory(nonempty_idx);
    tmp = tmp(nonempty_idx);
    nts = nts(nonempty_idx);
    meanmotion = meanmotion(nonempty_idx);

    include_idx = [];
    for i = 1:length(tmp)
        tt = tmp{i};
        if ~isempty(tt) & size(tt,1)>=375
            include_idx = [include_idx, i];
        end
    end

    tmp = tmp(include_idx);
    nts = nts(include_idx);
    meanmotion = meanmotion(include_idx);
    VisitIDvec_inventory = VisitIDvec_inventory(include_idx);
    sprintf('%d %d %d %d', size(tmp,1), size(nts,1), size(meanmotion,1), size(VisitIDvec_inventory,1))

    TSmat_pass = {};
    for idx = 1:length(tmp)
        TSmat_pass{idx} = tmp{idx}';
    end

    outputdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI';
    if exist(outputdir)
        mkdir(outputdir)
    end

    ROImat = struct;
    ROImat.TSmat_pass = TSmat_pass;
    ROImat.subjID = VisitIDvec_inventory;
    ROImat.nts = nts;
    ROImat.mean_motion = meanmotion;
    save(fullfile(outputdir, sprintf('%s_TSmat_pass.mat', imgname{taski})), 'ROImat', '-v7.3');
end

%%%%%% matching IDs
outputdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI';
rsfc_id = load(fullfile(outputdir, 'subjinfo.mat'));
rsfc_id = rsfc_id.subjinfo;

load(fullfile('/space/gwas-syn2/1/data/GWAS/ABCD/Images/NDA3.0/taskrsfc_R4/', sprintf('Release4.0_taskresid_Rs_roi_subinfo_baseline.mat')));
visitid_eventname = strcat(subjinfo.sitevec_inventory, '_', subjinfo.SubjIDvec_inventory, '_', subjinfo.eventname_inventory);
visitid_date = strcat(subjinfo.sitevec_inventory, '_', subjinfo.SubjIDvec_inventory, '_', subjinfo.datevec_inventory);

[c ia ib] = intersect(rsfc_id.subjID, visitid_date, 'stable');
visitid_eventname = visitid_eventname(ib);

subjinfo = rsfc_id;
subjinfo.mri_info_visitid = visitid_eventname;
save(fullfile(outputdir, 'subjinfo.mat'), 'subjinfo');
