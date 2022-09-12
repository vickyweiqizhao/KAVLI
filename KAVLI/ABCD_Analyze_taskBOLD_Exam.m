function ABCD_Analyze_taskBOLD_Exam(ContainerPath,FSContainerPath,varargin)
%function ABCD_Analyze_taskBOLD_Exam(ContainerPath,FSContainerPath,[options])
%
% Usage:
%  ABCD_Analyze_taskBOLD_Exam(ContainerPath,FSContainerPath,'key1', value1,...)
%
% Required Parameters:
%   ContainerPath: full path of BOLDPROC directory containing processed BOLD data
%   FSContainerPath: full path of directory containing freesurfer recon
%
% Mandatory Parameters:
%   'eprime_info': eprime info struct array
%     contains info about complete tfMRI scans
%       for all tasks in one scanning event
%     NOTE: must be supplied if fname_eprime is empty
%           will be used to set snums, runs, and fname_eprime
%     {default = []}
%   'full_info': protocol compliance & quality control info struct array
%     contains info about all scans in a scanning event
%     NOTE: must be supplied along with eprime_info
%     {default = []}
%   'fname_eprime': full path of csv or txt file derived from eprime output
%     may be cell array with number of elements matching snums
%     NOTE: must be supplied if eprime_info and full_info are empty
%     {default = []}
%   'snums': vector of scan numbers on which to run GLM analysis
%     NOTE: must be supplied if eprime_info and full_info are empty
%     {default = []}
%   'runs': vector of eprime run numbers corresponding to snums
%     if empty, will expect number of scans to match number of runs in eprime
%     {default = []}
%
% Optional Parameters:
%   'infix': BOLD file name infix (e.g. '', 'corr', 'corr_resBOLD')
%     {default = []}
%   'outstem': output file stem used for output directories
%     if empty, will construct from 'BOLD' and taskname
%     {default = []}
%   'analysis_outfix': string attached to analysis output subdirectories
%     {default = 'analysis'}
%   'max_eprime_delay_range': maximum eprime delay range (minutes)
%     {default = 1}
%   'max_eprime_delay': maximum eprime delay (minutes)
%     {default = 12.5}
%   'max_session_delay': maximum session delay (minutes)
%     for behavioral-only runs and for scans with consistent timing across runs
%     {default = 720}
%   'session_timing_flag': [0|1] allow delays over max
%     if there is consistent timing across runs and tasks
%     {default = 1}
%   'irregular_flag': [0|1|2] whether and how to filter output
%      0: all pGUID-event-tasks
%      1: pGUID-event-tasks with irregular acquisitions
%      2: pGUID-event-tasks with valid acquisition and timing
%     {default = 2}
%   'standard_flag': [0|1] whether to require standard acquisition
%     to not be considered "irregular" (2 scans, 1 eprime file) 
%     0: do not require standard acquisition
%     1: require standard acquisition after selecting matching scans
%     2: require standard acquisition before and after selecting matching scans
%     {default = 0}
%   'pguid_match_flag': [0|1] whether to consider irregular
%      if eprime_info.eprime_pguid_match = 0 for selected series
%     {default = 1}
%   'behav_only_flag': [0|1] whether to consider irregular
%      if eprime_info.eprime_behav_only = 1 for selected series
%     {default = 1}
%   'stim_times_flag': [0|1|2|3] use stimulus time files
%     0 = .1D (old style)
%     1 = .txt (model each TR of trial as an event)
%     2 = _block.txt and _dur.txt (model onset and duration of each event)
%     3 = _block.txt and dur = 0 (model onset of each event with impulse)
%     {default = 3}
%   'stim_times_model': name of stimulus model
%     model names in 3dDeconvolve for -stim_times
%     allowed: 'SPMG','GAM','TENT', 'BLOCK'
%     {default = 'SPMG'}
%   'stim_times_nbasis': number of basis functions for TENT
%     {default = 10}
%   'reml_flag': [0|1] whether to use AFNI's 3dREMLfit for pre-whitening
%     {default = false}
%   'concat_flag': [0|1|2] analyze concatenated across scans
%     0: analyze each scan individually
%     1: analyze concatenated scans only
%     2: analyze individually and concatenated
%     {default = 1}
%   'skipTRs': number of TRs at beginning of each run to ignore
%     {default = 0}
%   'minfrac': minimum fraction of a TR to register an event
%     {default = 0.5}
%   'minlag': number of TRs for minimum "lag" between stimulus and response
%     {default = 0}
%   'maxlag': number of TRs for maximum "lag" between stimulus and response
%     {default = 14}
%   'pthresh': probability threshold applied to f-stats for each contrast
%     if set to 0 or 1, no thresholding applied
%     {default = 0.05}
%   'taskname': name of fMRI task (e.g. MID, SST, nBack)
%     if empty, will construct from eprime filename
%     {default = []}
%   'contrasts_flag': [0|1] calculate glt contrasts between each condition
%     {default = 0}
%   'iresp_flag': [0|1] output impulse response functions for each condition
%     {default = 0}
%   'ico_flag': resample data to icosahedral sphere
%     {default = 1}
%   'ico_order': icosahedral order (0-7)
%      Order  Number of Vertices
%        0              12
%        1              42
%        2             162
%        3             642
%        4            2562
%        5           10242
%        6           40962
%        7          163842
%     {default = 5}
%   'ico_trunc_flag': [0|1] for ico_order<7, truncate extra values
%      otherwise, actually resample to ico_order
%     {default = 0}
%   'ico_presmooth': number of smoothing steps applied on native surface
%     before resampling to ico
%     NOTE: FWHM ~ 1.25 * sqrt(N)
%     {default = 16}
%   'ico_postsmooth': number of smoothing steps applied on ico surface
%     NOTE: FWHM ~ 1.25 * sqrt(N)
%     {default = 0}
%   'errts_flag': [0|1] whether to output the residuals in the deconvolution 
%     {default = true}
%   'vox_errts_flag': [0|1] whether to output the residuals in the deconvolution
%     in the voxel-native analysis 
%     {default = false}
%   'snums_valid': vector of scan numbers that were processed
%     if empty, will assume all were processed
%     {default = []}
%   'mc_flag': [0|1] whether within-scan motion correction was done
%     Will use motion.1D files as nuisance regressors for GLM
%     {default = 1}
%   'mc_inter_flag': [0|1] whether between-scan motion correction was done
%     Allows for a single reference scan to be used for registration to FS
%     {default = 1}
%   'regFS_flag': [0|1] whether to register BOLD scans to FreeSurfer nu.mgz
%     if registration already done, will not redo
%     necessary for painting to surface
%     if 0 and not already done, paint_flag is ignored
%     {default = 0}
%   'vox_flag': [0|1] perform voxel-wise GLM
%     {default = 1}
%   'mni_flag': [0|1] output voxel-wise data warped to MNI space
%     {default = 0}
%   'surf_flag': [0|1] perform surface vertex-wise GLM
%     {default = 0}
%   'roi_flag': [0|1] perform ROI-average GLM
%     {default = 0}
%   'outliers_flag': [0|1] censor outliers based on sem rms
%     {default = 0}
%   'outliers_thresh': sem rms threshold for censoring outliers  
%     {default = 5}
%   'aseg_flag': [0|1] extract ROI averages for subcortical segmentation
%     {default = 0}
%   'aseg_erode_flag': [0|1] "erode" aseg ROIs by smoothing and thresholding
%     after resampling to BOLD resolution
%     {default = 0}
%   'aseg_erode_nvoxels': smoothing kernel sigma for aseg erosion (# of voxels)
%     {default = 1}
%   'fname_aseg': name of aseg file
%      if empty, will use aseg.mgz in subjdir/subj/mri
%     {default = []}
%   'aparc_flag': [0|1] extract ROI averages for cortical surface parcellation
%     {default = 0}
%   'aparc_infix': string in aparc annot file
%      may be cell array of multiple infix strings
%      e.g. 'aparc', 'aparc.a2009s'
%     {default = 'aparc'}
%   'fparc_flag': extract time series for fparc cortical surface ROIs
%     ignored if fnames_fparc and fname_points are empty
%     {default = 1}
%   'fnames_fparc': cell array of annotation files in fsaverage space
%     will be resampled to to individual subject space before use
%     {default = []}
%   'render_flag': [0|1] whether to generate images of surface maps
%     {default = 0}
%   'projdist': dist (mm) to project along normal when painting
%     {default = 1}
%   'projfrac': fractional dist to project along normal when painting
%     {default = 0.5}
%   'projfrac_flag': [0|1] whether to use projfrac (1) or projdist (0)
%     {default = 0}
%   'projfrac_avg': vector of [min max del] for averaging multiple samples
%     If empty, use projfrac instead if projfrac_flag=1
%     {default = []}
%   'projdist_avg': vector of [min max del] for averaging multiple samples
%     If empty, use projdist instead
%     {default = []}
%   'resamp_flag': [0|1] whether to resample results to 1x1x1 before painting
%     {default = 0}
%   'force_repaint_flag': [0|1] whether to repaint even if output files exist
%     (do not redo volume calculations)
%     {default = 0}
%   'surfname': name of surface onto which to sample volume data
%     {default = white}
%   'verbose': [0|1|2] display status messages
%      0: no messages except errors
%      1: no messages except WARNING
%      2: frequent status messages
%     {default = 1}
%   'check_stale_flag': check creation date of output directory
%     and delete if stale (created prior to fs.finish.all.touch)
%     {default = 1}
%   'check_complete_flag': [0|1] whether to require that recon is complete
%     {default = 1}
%   'FS_version': which version of Freesurfer used (e.g. 305, 450, 510, 530)
%     for checking whether recon is complete and loading FreeSurferColorlUT
%     {default = 530}
%   'stepstr': string used to record completed processing
%     {default = 'analyzing_taskbold'}
%   'forceflag': [0|1] whether to run calculations even if output files exist
%     {default = 0}
%
% Created:  12/02/16 by Don Hagler
% Prev Mod: 08/06/19 by Dani Cornejo
% Prev Mod: 04/01/21 by Don Hagler
% Last Mod: 06/17/21 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not modify this section:

% based on MMIL_Analyze_taskBOLD_Exam
% Created:  09/01/08 by Don Hagler
% Last Mod: 03/07/13 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;

% check input parameters
[parms,errcode] = check_input(ContainerPath,FSContainerPath,varargin);
if errcode, write_progress(parms,errcode); return; end;

% get number of scans, etc.
[parms,errcode] = check_sess(parms);
if errcode, write_progress(parms,errcode); return; end;

% check freesurfer, regFS
[parms,errcode] = check_fsurf(parms);
if errcode, write_progress(parms,errcode); return; end;

% get file names for each scan
[parms,errcode] = check_scans(parms);
if errcode, write_progress(parms,errcode); return; end;

% create stimulus time course files
[parms,errcode] = prep_stims(parms); 
if errcode, write_progress(parms,errcode); return; end;

% create condition contrast structure
parms = prep_contrasts(parms);

% set stimulus file names
parms = set_stim_fnames(parms);

% create cell array with sets of scan numbers
parms = set_snum_sets(parms);

% create output directory / directories
parms = prep_outdir(parms);

% resample fparc ROIs to subject
parms = resamp_fparc(parms);

% run GLM analysis
[parms,errcode] = run_analysis(parms);
if errcode, write_progress(parms,errcode); return; end;

% write progress file
write_progress(parms,errcode);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_progress(parms,errcode)
  mmil_write_progress(parms.cpath,parms.step,errcode);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [parms,errcode] = check_input(ContainerPath,FSContainerPath,options)
  errcode = 0;
  parms = mmil_args2parms(options,{...
    'cpath',ContainerPath,[],...
    'fspath',FSContainerPath,[],...
  ...
    'fname_eprime',[],[],...
    'snums',[],[],...
    'runs',[],[],...
    'eprime_info',[],[],...
    'full_info',[],[],...
  ...
    'infix',[],[],...
    'outstem',[],[],...
    'analysis_outfix','analysis',[],...
    ...
    'max_eprime_delay_range',1,[],...
    'max_eprime_delay',12.5,[],...
    'max_session_delay',720,[],...
    'session_timing_flag',true,[false true],...
    'ignore_delay_flag',false,[false true],...
    'irregular_flag',2,[0,1,2],...
    'standard_flag',0,[0,1,2],...
    'pguid_match_flag',true,[false true],...
    'behav_only_flag',true,[false true],...
    ...
    'stim_times_flag',3,0:3,...
    'stim_times_model','SPMG',{'SPMG','TENT','GAM','BLOCK'},...
    'stim_times_nbasis',10,[1,100],...
    'reml_flag',false,[false true],...
    'concat_flag',1,0:2,...
    'contrasts_flag',false,[false true],...
    'iresp_flag',false,[false true],...
    'skipTRs',0,[0,1000],...
    'minfrac',0.5,[0,1],...
    'minlag',0,[0,10],...
    'maxlag',14,[0,30],...
    'tksmooth',10,[0,1000],...
    'pthresh',0.05,[0,1],...
    'taskname',[],[],...
    'ico_flag',true,[false true],...
    'ico_order',5,[0 7],...
    'ico_trunc_flag',false,[false true],...
    'ico_presmooth',16,[0,1000],...
    'ico_postsmooth',0,[0,1000],...
    'stats',{'beta','sem'},{'beta','sem'},...
    'errts_flag',true,[false true],...
    'vox_errts_flag',false,[false true],...
  ...
    'stim_fnames',[],[],...
    'stim_labels',[],[],...
    'snums_valid',[],[],...
    'mc_flag',true,[false true],...
    'mc_inter_flag',true,[false true],...
    'regFS_flag',false,[false true],...
    'vox_flag',true,[false true],...
    'mni_flag',false,[false true],...
    'surf_flag',false,[false true],...
    'roi_flag',false,[false true],...
    'outliers_flag',false,[false true],...
    'outliers_thresh',5,[0,10],...
    'paint_flag',true,[false true],...
    'aseg_flag',false,[false true],...
    'aseg_erode_flag',false,[false true],...
    'aseg_erode_nvoxels',1,[1:100],...
    'fname_aseg',[],[],...
    'aparc_flag',false,[false true],...
    'aparc_infix','aparc',[],...
    'fparc_flag',true,[false true],...
    'fnames_fparc',[],[],...
    'render_flag',false,[false true],...
    'projdist',1,[0,10],...
    'projfrac',0.5,[0,2],...
    'projfrac_flag',false,[false true],...
    'projdist_avg',[],[],...
    'projfrac_avg',[],[],...
    'paint_surf','white',[],...
    'resamp_flag',false,[false true],...
    'force_repaint_flag',false,[false true],...
    'surfname','white',[],...
    'mask_midbrain_flag',false,[false true],...
    'verbose',1,[0:2],...
    'check_stale_flag',true,[false true],...
    'check_complete_flag',true,[false true],...
    'FS_version',530,[],...
    'stepstr','analyzing_taskbold',[],...
    'forceflag',false,[false true],...    
  ... % parameters for rendering output
    'curvflag',true,[false true],...
    'curvfact',0.2,[0,1],...
    'zoom',1,[0.01,100],...
    'tif_flag',true,[false true],...
    'tif_dpi',300,[10,1000],...
    'fmin',0.1,[0,100],...
    'fmid',1,[0.001,1000],...
    'fmax',2,[0.001,1000],...
    'view_surfname','inflated',{'inflated','white','pial'},...
    'viewlist',{'lat','med'},{'left','right','lat','med','pos',...
                              'ant','sup','ven','infr','infl'},...
  ...
    'res_outfix','resBOLD',[],...
  ...
    'motion_radius',50,[],... % for calculating distance from angle
    'motion_absflag',true,[false true],...
    'motion_nodyflag',false,[false true],...
    'censor_flag',true,[false true],... 
    'censor_thresh',0.9,[0.001,Inf],...
    'contig_thresh',0,[0,20],...
    'mc_resp_flag',true,[false true],...
    'resp_low',18.6,[],...
    'resp_high',25.7,[],...
  ...
    'min_numTRs',300,[],...
    'fnamestem','BOLD',[],...
    'out_ext','.mgz',{'.mgh','.mgz'},...
    'revflag',false,[false true],... % only applys for pep and ipp
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'ico_nverts',[42 162 642 2562 10242 40962 163842],[],...
  ...
    'skipTRs_GE',5,[],...
    'skipTRs_GE_DV26',16,[],...
    'skipTRs_Philips',8,[],...
    'skipTRs_Siemens',8,[],...
  ... % overridden later
    'conds_contrast',[],[],...
  ... % parameter names
    'paint_tags',{'outtype','regfile' 'projfrac_flag' 'projdist' 'projfrac'...
                  'projdist_avg' 'projfrac_avg' 'mask_midbrain_flag' ...
                  'subjdir' 'surfname' 'overwrite_flag'},[],...
    'sview_tags',{'frames','subjdir','surfname','surfdir','hemi',...
                  'colorscale','cmap','fmax','fmin','fmid',...
                  'curvflag','curvfact','view','zoom',...
                  'tif_flag','tif_dpi','outstem','outdir',...
                  'visible_flag','pause_dur','forceflag'},[],...
    'deconv_tags',{'stim_labels','stim_times_flag','stim_times_model',...
                   'stim_times_nbasis','contrasts_flag','iresp_flag',...
                   'censor_flag','censor_thresh','contig_thresh',...
                   'deriv_flag','reml_flag','errts_flag','vox_errts_flag',...
                   'outdir','outstem','skipTRs',...
                   'TR','norm_flag','detrend','thresh','out_ext','nii_ext',...
                   'minlag','maxlag','fname_motion','motion_labels',...
                   'glt_fnames','glt_labels','cleanupflag','forceflag',...
                   'goforit_flag','goforit_val','vox_flag','surf_flag','roi_flag',...
                   'motion_radius','motion_absflag','motion_nodyflag',...
                   'mc_resp_flag','resp_low','resp_high','conds_contrast',...
                   'regfile','subjdir','subj','ico_presmooth','ico_flag',...
                   'ico_order','ico_trunc_flag','ico_postsmooth',...
                   'aseg_flag','aseg_erode_flag','aseg_erode_nvoxels',...
                   'fname_aseg','aparc_flag','fnames_aparc','num_aparcs',...
                   'aparc_names','aparc_hemis','hemilist'},[],...  
    'aseg_tags',{'outdir','outstem','fname_out','csv_flag','fname_aseg',...
                 'aseg_aparc_flag','fname_vals','dispvec',...
                 'disp_roicodes','disp_scalefact','dispfact',...
                 'erode_flag','erode_nvoxels',...
                 'scalefact','minval','M_reg','res_outfix',...
                 'FS_version','fname_colorlut',...
                 'verbose','forceflag',...
                 'aseg_roilist','aparc_roilist','exclude_roilist','frames'},[],...
    'surf_roi_tags',{'fname_aparc','fname_label','fname_weights','frames',...
                     'minval','scalefact','fname_colorlut','hemi',...
                     'weights_thresh','verbose','annot_name'},[],...
    'skip_tags',{'skipTRs_GE','skipTRs_GE_DV26',...
                'skipTRs_Philips','skipTRs_Siemens'},[],...
    'runs_tags',{'max_eprime_delay_range','max_eprime_delay','max_session_delay',...
                'ignore_delay_flag','session_timing_flag','irregular_flag','standard_flag',...
                'pguid_match_flag','behav_only_flag','verbose'},[],...
    'outdir_tags',{'runs','outstem','fnamestem',...
                   'infix','analysis_outfix'},[],...
  });
  
  parms.nhemi = length(parms.hemilist);
  % check input directories
  if isempty(parms.cpath)
    error('ContainerPath is empty');
  elseif ~exist(parms.cpath,'file')
    error('ContainerPath %s not found',parms.cpath);
  end;
  if isempty(parms.fspath)
    error('FSContainerPath is empty');
  elseif ~exist(parms.fspath,'file')
    error('FSContainerPath %s not found',parms.fspath);
  end;
  
  % check that fname_eprime supplied or eprime_info and full_info
  if isempty(parms.fname_eprime) &&...
     (isempty(parms.eprime_info) || isempty(parms.full_info))
    fprintf('%s: ERROR: no E-prime info supplied\n',mfilename);
    errcode = 1;
    return;
  end;

  % get task name if not supplied
  if isempty(parms.taskname)
    if ~isempty(parms.eprime_info)
      % get task name from first entry in eprime_info
      parms.taskname = parms.eprime_info(1).taskname;
    else
      % set task name based on eprime file name
      parms.taskname = get_taskname(parms.fname_eprime);
    end;
  end

  % parameters for fs_paint
  parms.overwrite_flag = (parms.forceflag || parms.force_repaint_flag);
  % mmil_aseg_analysis
  parms.erode_flag = parms.aseg_erode_flag;
  parms.erode_nvoxels = parms.aseg_erode_nvoxels;
  % set processing step name
  parms.step = sprintf('%s_%s',parms.stepstr,lower(parms.taskname));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function taskname = get_taskname(fname_eprime)
  taskname = [];
  if iscell(fname_eprime), fname_eprime = fname_eprime{1}; end;
  [fpath,fstem] = fileparts(fname_eprime);
  if ~isempty(regexpi(fstem,'MID')) 
    taskname = 'MID';
  elseif ~isempty(regexpi(fstem,'SST')) 
    taskname = 'SST';
  elseif ~isempty(regexpi(fstem,'WM')) 
    taskname = 'nBack';
  else
    k = regexp(fstem,'\w+_(?<task>\w+)$','names');
    taskname = k.task;
  end;
  if strcmp(upper(taskname),'WM')
    taskname = 'nBack';
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [parms,errcode] = check_sess(parms)
  errcode = 0;
  
  % load ContainerInfo
  fname_info = sprintf('%s/ContainerInfo.mat',parms.cpath);
  if ~exist(fname_info,'file')
    fprintf('%s: ERROR: missing ContainerInfo file %s\n',...
        mfilename,fname_info);
    errcode = 1;
    return;
  end;
  [ContainerInfo,errcode] = MMIL_Load_ContainerInfo(parms.cpath);
  if errcode, return; end;

  % check task fMRI scans and runs
  if isempty(parms.fname_eprime) &&...
    ~isempty(parms.eprime_info) && ~isempty(parms.full_info)
    if parms.verbose>1
      fprintf('%s: checking tfMRI runs...\n',mfilename);
    end;
    tmp_parms = parms;
    tmp_parms.verbose = (parms.verbose>1);
    a = mmil_parms2args(tmp_parms,parms.runs_tags);
    try
      tfMRI_info = abcd_check_tfMRI_runs(parms.eprime_info,parms.full_info,ContainerInfo,a{:});
    catch me
      fprintf('%s: ERROR: abcd_check_tfMRI_runs failed:\n%s\n',...
        mfilename,me.message);
      errcode = 1;
      return;
    end;
    % get tfMRI_info for this task
    idx_task = find(strcmp(parms.taskname,{tfMRI_info.taskname}));
    tfMRI_info = tfMRI_info(idx_task);
    if isempty(tfMRI_info)
      fprintf('%s: ERROR: no usable %s data\n',...
        mfilename,parms.taskname);
      errcode = 1;
      return;
    end;
    % set parameters for scans selected and found
    parms.snums = tfMRI_info.snums;
    parms.runs = tfMRI_info.runs;
    parms.eprime_runs = tfMRI_info.eprime_runs;
    if parms.verbose>1
      fprintf('%s: snums = [ %s], runs = [ %s], eprime_runs = [ %s]\n',...
        mfilename,sprintf('%d ',parms.snums),sprintf('%d ',parms.runs),...
                  sprintf('%d ',parms.eprime_runs));
    end;
    parms.fname_eprime = tfMRI_info.fnames_eprime;
  else
    if ~iscell(parms.fname_eprime)
      parms.fname_eprime = {parms.fname_eprime};
    end;
    if length(parms.fname_eprime)==1 && ...
       length(parms.fname_eprime)~=length(parms.snums)
      parms.fname_eprime = repmat(parms.fname_eprime,[1,length(parms.snums)]);
    end;
    if length(parms.fname_eprime)~=length(parms.snums)
      fprintf('%s: ERROR: mismatch between fname_eprime and snums\n',...
        mfilename);
      errcode = 1;
      return;
    end;
  end;
  if isempty(parms.fname_eprime)
    fprintf('%s: ERROR: empty fname_eprime\n',mfilename);
    errcode = 1;
    return;
  end;
  if parms.verbose>1
    fprintf('%s: fname_eprime = \n%s',...
      mfilename,sprintf('%s\n',parms.fname_eprime{:}));
  end;
  
  % set number of TRs to skip
  %   based on Manufacturer and SoftwareVersion
  args = mmil_parms2args(parms,parms.skip_tags);
  parms.skipTRs = abcd_set_skipTRs(ContainerInfo,args{:});
  if parms.verbose>1
    fprintf('%s: NOTE: skipTRs = %d\n',mfilename,parms.skipTRs);
  end;

  % get number of scans, etc.
  [parms.ScanInfo,parms.SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(parms.cpath,...
    'snums',parms.snums_valid,'fnamestem',parms.fnamestem);
  if errcode || isempty(parms.ScanInfo)
    errcode = 1;
    return;
  end;
  fprintf('%s: %d BOLD scans in %s\n',...
    mfilename,parms.SessInfo.nscans,parms.cpath);
  if isempty(parms.snums)
    parms.snums = [1:parms.SessInfo.nscans];
  end;
  if isempty(parms.snums)
    fprintf('%s: ERROR: no valid scans\n',mfilename);
    errcode = 1;
    return;
  end;

  % check runs
  if ~isempty(parms.runs) && length(parms.runs) ~= length(parms.snums)
    fprintf('%s: ERROR: length of runs (%d) does not match length of snums (%d)\n',...
      mfilename,length(parms.runs),length(parms.snums));
    errcode = 1;
    return;
  end;
  % check eprime_runs
  if ~isempty(parms.eprime_runs) && length(parms.eprime_runs) ~= length(parms.snums)
    fprintf('%s: ERROR: length of eprime_runs (%d) does not match length of snums (%d)\n',...
      mfilename,length(parms.eprime_runs),length(parms.snums));
    errcode = 1;
    return;
  end;
  
  % check for more than 2 runs
  if length(parms.snums)>2
    fprintf('%s: WARNING: limiting to 2 runs (instead of %d)\n',...
      mfilename,length(parms.snums));
    parms.snums = parms.snums(1:2);
    if ~isempty(parms.runs)
      parms.runs = parms.runs(1:2);
    end;
    if ~isempty(parms.eprime_runs)
      parms.eprime_runs = parms.eprime_runs(1:2);
    end;
  end;

  % check for duplicate run numbers
  if length(parms.runs)>1 && length(unique(parms.runs))<length(parms.runs)
    if strcmp(parms.taskname,'SST')
      % runs are interchangeable for SST, so number them sequentially
      fprintf('%s: WARNING: treating duplicate run numbers as sequential because task is %s\n',...
        mfilename,parms.taskname);
      parms.runs = [1:length(parms.snums)];
    else
      % runs are NOT interchangeable for MID and nBack, so exclude 2nd run
      if all(parms.runs==1)
        fprintf('%s: WARNING: limiting to first run because of multiple run 1 with task %s\n',...
          mfilename,parms.taskname);
        parms.snums = parms.snums(1);
        parms.runs = parms.runs(1);
        if ~isempty(parms.eprime_runs)
          parms.eprime_runs = parms.eprime_runs(1);
        end;
      else
        fprintf('%s: ERROR: duplicate run numbers\n',mfilename);
        errcode = 1;
        return;
      end;
    end;
  end;
  
  % get nreps and TRs
  parms.nreps = [parms.ScanInfo(parms.snums).nreps];   
  parms.TRs = [parms.ScanInfo(parms.snums).TR]/1000;
  % exclude scans with nreps < min_numTRs, but use the rest
  ind_keep = find(parms.nreps>parms.min_numTRs); 
  if isempty(ind_keep)
    fprintf('%s: ERROR: numTRs below minimum (%d): %d: scans %s\n',...
      mfilename,parms.min_numTRs,parms.nreps,strtrim(sprintf('%d ',parms.snums)));
    errcode = 1; return;
  end;
  parms.snums = parms.snums(ind_keep);
  if ~isempty(parms.runs)
    parms.runs = parms.runs(ind_keep);
  end;
  if ~isempty(parms.eprime_runs)
    parms.eprime_runs = parms.eprime_runs(ind_keep);
  end;
  parms.nreps = parms.nreps(ind_keep);
  parms.TRs = parms.TRs(ind_keep);
  % use maximum numTRs and TR if multiple scans 
  if length(parms.snums)>1
    % check for mismatch in numTRs
    if length(unique(parms.nreps))>1
      fprintf('%s: WARNING: mismatch in numTRs:\n',mfilename);
      for i=1:length(parms.snums)
        fprintf('  scan %d: %d TRs\n',...
          parms.snums(i),parms.nreps(i));
      end;
    end;
    % check for mismatch in TRs
    if length(unique(parms.TRs))>1
      fprintf('%s: WARNING: mismatch in TRs:\n',mfilename);
      for i=1:length(parms.snums)
        fprintf('  scan %d: %d TRs\n',...
          parms.snums(i),parms.TRs(i));
      end;
    end;
  end;
  % use maximum numTRs and TR if multiple scans
  parms.numTRs = max(parms.nreps);
  parms.TR = max(parms.TRs);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_fsurf(parms)
  errcode = 0;
  % for ROI analysis
  if parms.aparc_flag || parms.fparc_flag, parms.paint_flag = 1; end;
  % check FreeSurfer path and regT1
  if parms.aseg_flag || parms.paint_flag || parms.resamp_flag
    if isempty(parms.fspath)
      fprintf('%s: WARNING: setting aseg_flag and paint_flag = 0 because fspath is empty',mfilename);
    elseif ~exist(parms.fspath,'dir')
      fprintf('%s: WARNING: setting aseg_flag and paint_flag = 0 because fspath %s not found',...
        mfilename,parms.fspath);
    end;
    if isempty(parms.fspath) || ~exist(parms.fspath,'dir')
      parms.resamp_flag = 0; parms.aseg_flag = 0;
      parms.paint_flag = 0; parms.aparc_flag = 0; parms.fparc_flag = 0;
      return;
    end;
    % registration to T1
    if parms.regFS_flag
      fname_T1 = sprintf('%s/mri/nu.mgz',parms.fspath);
      BOLD_MMIL_Register_to_T1(...
        parms.cpath,...
        'fname_T1',fname_T1,...
        'snums',parms.snums_valid,...
        'infix',parms.infix,...
        'forceflag',parms.forceflag);
    end;
    % get subjdir from fspath
    [parms.subjdir,parms.subj,ext] = fileparts(parms.fspath);
    parms.subj = [parms.subj ext]; % in case of .'s and such
    % check status of FreeSurfer recon
    if parms.check_complete_flag
      [status,message] = MMIL_Get_FSReconStatus(parms.fspath,parms.FS_version);
      if ~ismember(status,[2,5,6])
        fprintf('%s: WARNING: incomplete recon for %s\n',mfilename,parms.subj);
        errcode = 1;
        return;
      end;
      if status==6
        fprintf('%s: WARNING: only volume recon is complete for %s\n',...
          mfilename,parms.subj);
        parms.paint_flag = 0;
        parms.aparc_flag = 0;
        parms.fparc_flag = 0;
      end;
    end;
    setenv('SUBJECTS_DIR',parms.subjdir); % NOTE: not sure this is necessary
    fprintf('%s: using scan %d as reference\n',...
      mfilename,parms.SessInfo.regT1_ref);
    parms.fstem_ref = [parms.cpath '/' parms.SessInfo.fstem_regT1_ref];
    if ~isempty(parms.infix)
      parms.fstem_ref = [parms.fstem_ref '_' parms.infix];
    end;
  end;
  % check fname_aseg
  if parms.aseg_flag
    if isempty(parms.fname_aseg)
      parms.fname_aseg = sprintf('%s/mri/aseg.mgz',parms.fspath);
    end;
    if ~exist(parms.fname_aseg,'file')
      fprintf('%s: ERROR: aseg file %s not found\n',mfilename,parms.fname_aseg);
      errcode = 1;
      return;
    end;
  else
    parms.fname_aseg = [];
  end;
  % check whether aparc_infix is a cell array
  if parms.aparc_flag && ~iscell(parms.aparc_infix)
    parms.aparc_infix = {parms.aparc_infix};
  end;
  % set fnames_aparc
  parms.aparc_hemis = [];
  parms.aparc_names = [];
  parms.fnames_aparc = [];
  parms.num_aparcs = 0;
  if parms.aparc_flag
    parms.naparc = length(parms.aparc_infix);
    k = 1;
    for i=1:parms.naparc
      aparc_infix = parms.aparc_infix{i};
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        parms.aparc_hemis{k} = hemi;
        parms.aparc_names{k} = aparc_infix;
        parms.fnames_aparc{k} = sprintf('%s/label/%s.%s.annot',...
          parms.fspath,hemi,aparc_infix);
        k = k + 1;
      end;
    end;
    parms.num_aparcs = length(parms.fnames_aparc);
    for f=1:parms.num_aparcs
      if ~exist(parms.fnames_aparc{f},'file')
        fprintf('%s: ERROR: aparc annot file %s not found\n',...
          mfilename,parms.fnames_aparc{f});
        errcode = 1;
        return;
      end;
    end;
  end;
  % check fnames_fparc
  if parms.fparc_flag && ~isempty(parms.fnames_fparc)
    if ~iscell(parms.fnames_fparc)
      parms.fnames_fparc = {parms.fnames_fparc};
    end;
    if length(parms.fnames_fparc) ~= length(parms.hemilist)
      error('number of elements in fnames_fparc must match hemilist');
    end;
    for f=1:length(parms.fnames_fparc)
      if ~exist(parms.fnames_fparc{f},'file')
        fprintf('%s: ERROR: fparc annot file %s not found\n',...
          mfilename,parms.fnames_fparc{f});
        errcode = 1;
        return;
      end;
    end;    
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_aparc = set_fname_aparc(parms,hemi,aparc_infix)
  if strcmp(aparc_infix,'fparc')
    h = find(strcmp(hemi,parms.fparc_hemis));
    fname_aparc = parms.fnames_fparc_resamp{h};
  else
    fname_aparc = sprintf('%s/label/%s.%s.annot',...
      parms.fspath,hemi,aparc_infix);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = check_scans(parms)
  errcode = 0;
  % check input, save input file names
  parms.data_fnames = [];
  parms.data_fstems = [];
  parms.reg_fnames = [];
  parms.motion_fnames = [];
  for i=1:length(parms.snums)
    s = parms.snums(i);
    % check for bad scan num
    if ~ismember(s,parms.SessInfo.snums_valid)
      fprintf('%s: ERROR: bad BOLD Scan Num (%d)\n',mfilename,s);
      errcode = 1;
      return;
    end;
    % check motion and data files exist
    fstem = sprintf('%s/%s',parms.cpath,parms.ScanInfo(s).fstem);
    if ~isempty(parms.infix)
      fstem = [fstem '_' parms.infix];
    end;
    if parms.mc_flag
      parms.motion_fnames{i} = [fstem '_motion.1D'];
      if ~exist(parms.motion_fnames{i},'file')
        fprintf('%s: ERROR: motion 1D file %s not found\n',...
          mfilename,parms.motion_fnames{i});
        errcode = 1;
        return;
      end;
    end;
    parms.data_fstems{i} = fstem;
    parms.data_fnames{i} = sprintf('%s.mgz',fstem);
    if ~exist(parms.data_fnames{i},'file')
      fprintf('%s: ERROR: data file %s not found\n',mfilename,parms.data_fnames{i});
      errcode = 1;
      return;
    end;
    if parms.paint_flag || parms.resamp_flag
      % check register.dat exists
      if parms.mc_inter_flag
        parms.reg_fnames{i} = sprintf('%s_register.dat',parms.fstem_ref);
      else
        parms.reg_fnames{i} = sprintf('%s_register.dat',fstem);
      end;
      if ~exist(parms.reg_fnames{i},'file')
        fprintf('%s: WARNING: reg file %s not found (cannot resample to T1 or paint to surface)\n',...
          mfilename,parms.reg_fnames{i});
        parms.resamp_flag = 0; parms.aseg_flag = 0;
        parms.paint_flag = 0; parms.aparc_flag = 0;
        parms.fparc_flag = 0;
      end;
    end;
  end;
  parms.nscans = length(parms.snums);
  if ~parms.nscans
    fprintf('%s: ERROR: no scans found\n',mfilename);
    errcode = 1;
    return;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = prep_stims(parms)
  errcode = 0;
  % set parameters for extracting from eprime files
  tp = [];
  tp.outstem = [];
  tp.outdir = [parms.cpath '/stim1D_' parms.taskname];
  tp.numTRs = parms.numTRs;
  tp.TR = parms.TR;
  tp.minfrac = parms.minfrac;
  tp.nskipTRs = 0; % stimuli start after skipped TRs
  tp.forceflag = parms.forceflag;
  args = mmil_parms2args(tp);
  % create 1D files for each unique eprime file
  fname_eprime = unique(parms.fname_eprime);
  parms.eprime_nruns = 0;
  for i=1:length(fname_eprime)
    switch upper(parms.taskname)
      case 'MID'
        [eprime_nruns,errcode] = abcd_extract_eprime_mid(fname_eprime{i},args{:});
      case 'SST'
        [eprime_nruns,errcode] = abcd_extract_eprime_sst(fname_eprime{i},args{:});
      case 'NBACK'
        [eprime_nruns,errcode] = abcd_extract_eprime_nback(fname_eprime{i},args{:});
      otherwise
        error('prep_stims not implemented for %s',parms.taskname);
    end;
    parms.eprime_nruns = parms.eprime_nruns + eprime_nruns;    
  end;
  parms.stimdir = tp.outdir;
  if errcode, return; end;
  if parms.verbose>1
    fprintf('%s: exporting events...\n',mfilename);
  end;
  for s=1:parms.nscans
    fname_eprime = parms.fname_eprime{s};
    if ~isempty(parms.runs)
      run = parms.runs(s);
    else
      run = s;
    end;
    if ~isempty(parms.eprime_runs)
      eprime_run = parms.eprime_runs(s);
    else
      eprime_run = s;
    end;
    abcd_export_events(parms.cpath,parms.taskname,fname_eprime,run,eprime_run,parms.forceflag);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_stim_fnames(parms)
  if any(~ismember(parms.snums,1:parms.SessInfo.nscans))
    error('bad scan numbers in snums');
  end;
  nscans = length(parms.snums);
  fprintf('%s: number of scans = %d\n',mfilename,nscans);
  if nscans==0, return; end;
  fprintf('%s: run numbers = [%s]\n',mfilename,sprintf('%d ',parms.runs));
  fprintf('%s: eprime run numbers = [%s]\n',mfilename,sprintf('%d ',parms.eprime_runs));
  nstims = length(parms.stim_labels);
  fprintf('%s: number of stimulus conditions = %d\n',mfilename,nstims);
  parms.stim_fnames = cell(parms.nscans,1);
  switch parms.stim_times_flag
    case 0
      fext = '.1D';
    case 1
      fext = '.txt';
    case {2,3}
      fext = '_block.txt';
  end;
  % compile arrays of stimulus timing files for each run
  for s=1:parms.nscans
    [tmp,fstem_eprime] = fileparts(parms.fname_eprime{s});
    fstem_eprime = abcd_clean_fstem(fstem_eprime);
    stim_fnames = cell(nstims,1);
    if ~isempty(parms.eprime_runs)
      eprime_run = parms.eprime_runs(s);
    else
      eprime_run = s;
    end;
    for c=1:nstims
      stim_label = parms.stim_labels{c};
      stim_fname = sprintf('%s/%s_scan%d_%s%s',...
        parms.stimdir,fstem_eprime,eprime_run,stim_label,fext);
      if ~exist(stim_fname,'file')
        fprintf('%s: WARNING: stim file %s not found, creating null stim file...\n',... 
          mfilename,stim_fname);
        % create file with all zeros (1D) or * (txt)
        switch fext
          case '.txt'
            % write to txt file
            fid = fopen(stim_fname,'wt');
            if fid<0, error('failed to open %s for writing',stim_fname); end;
            fprintf(fid,'*\n');
            fclose(fid);
          case '.1D'
            % write to 1D file
            fid = fopen(stim_fname,'wt');
            if fid<0, error('failed to open %s for writing',stim_fname); end;
            for k=1:parms.numTRs
              fprintf(fid,'%d\n',0);
            end;
            fclose(fid);
        end;
      end;
      stim_fnames{c} = stim_fname;
    end;
    parms.stim_fnames{s} = stim_fnames;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_contrasts(parms)
  switch upper(parms.taskname)
    case 'MID'
      [parms.stim_labels,parms.conds_contrast] = abcd_set_contrasts_mid();
    case 'SST'
      [parms.stim_labels,parms.conds_contrast] = abcd_set_contrasts_sst();
    case 'NBACK'
      [parms.stim_labels,parms.conds_contrast] = abcd_set_contrasts_nback();
  end;
  parms.nstims = length(parms.stim_labels);
  parms.ncontrasts = length(parms.conds_contrast);
  parms.condnames = cat(2,parms.stim_labels,{parms.conds_contrast.name});
  parms.nconds = parms.nstims + parms.ncontrasts;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_snum_sets(parms)
  parms.snum_sets = {};
  if parms.nscans==1 || ismember(parms.concat_flag,[0,2])
    parms.snum_sets = num2cell(parms.snums);
  end;
  if parms.nscans>1 && parms.concat_flag>0
    parms.snum_sets{end+1} = parms.snums;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_outdir(parms)
  if isempty(parms.outstem)
    parms.outstem = sprintf('%s_%s',parms.fnamestem,parms.taskname);
  end;
  % create output dirs
  parms.outdirs = cell(size(parms.snum_sets));
  parms.outstems = cell(size(parms.snum_sets));
  for i=1:length(parms.snum_sets)
    s = parms.snum_sets{i};
    [~,idx_scan] = intersect(parms.snums,s);
    snums = parms.snums(idx_scan);
    tp = parms;
    tp.runs = parms.runs(idx_scan);
    args = mmil_parms2args(tp,parms.outdir_tags);
    [parms.outdirs{i},parms.outstems{i}] = ...
      abcd_taskBOLD_outdir(snums,parms.taskname,args{:});
    parms.outdirs{i} = [parms.cpath '/' parms.outdirs{i}];
    % check timestamp of final or surf touch file
    %   compare to timestamp of outdir
    check_stale(parms.outdirs{i},parms);
    mmil_mkdir(parms.outdirs{i});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_stale(outdir,parms)
  if parms.check_stale_flag && exist(outdir,'dir')
    fstem_touch = 'fs.finish.all';
    fname_touch = sprintf('%s/touch/%s.touch',parms.fspath,fstem_touch);
    if ~exist(fname_touch,'file')
      fprintf('%s: WARNING: %s not found for %s (unfinished recon)\n',...
        mfilename,fstem_touch,parms.subj);
      return;
    else
      dir_touch = dir(fname_touch);
      dir_outdir = dir(outdir);
      date_touch = dir_touch.datenum;
      date_outdir = min([dir_outdir.datenum]);
      if date_outdir < date_touch
        fprintf('%s: WARNING: fs.finish.all.touch is more recent than %s, removing existing output...\n',...
          mfilename,outdir);
        rmdir(outdir,'s');
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resample_fparc_annot(parms)
  i = parms.num_aparcs;
  for f=1:length(parms.fnames_fparc)
    fname_in = parms.fnames_fparc{f};
    if parms.verbose>1
      fprintf('%s: resampling annotation file %s from fsaverage to %s...\n',...
        mfilename,fname_in,parms.subj);
    end;
    % call fs_annot2annot
    tmp_parms = [];
    tmp_parms.outdir = [parms.cpath '/fparc_annot'];
    tmp_parms.source_subj = 'fsaverage';
    tmp_parms.subj = parms.subj;
    tmp_parms.subjdir = parms.subjdir;
    tmp_parms.verbose = (parms.verbose>1);
    tmp_parms.forceflag = parms.forceflag;
    args = mmil_parms2args(tmp_parms);
    % check timestamp of final or surf touch file, compare to timestamps within outdir
    check_stale(tmp_parms.outdir,parms);
    % resample annotation from ico to subj
    fname_out = fs_annot2annot(fname_in,args{:});
    n = regexp(fname_out,'(?<hemi>[lr]h)\.(?<name>.+)\.annot$','names');
    if isempty(n)
      error('unexpected naming convention for aparc file %s\n',fname_out);
    end;
    parms.fnames_fparc_resamp{f} = fname_out;
    parms.fparc_hemis{f} = n.hemi;
    parms.fparc_names{f} = n.name;
    i = i + 1;
    parms.fnames_aparc{i} = fname_out;
    parms.aparc_hemis{i} = n.hemi;
    parms.aparc_names{i} = n.name;
  end;
  parms.num_aparcs = length(parms.fnames_aparc);
  parms.aparc_flag = 1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = run_analysis(parms)
  errcode = 0;
  for i=1:length(parms.snum_sets)
    s = parms.snum_sets{i};
    fprintf('%s: running GLM analysis for',mfilename);
    if length(s)==1
      fprintf(' scan %d\n',s);
    else
      fprintf(' scans [%s]\n',strtrim(sprintf('%d ',s)));
    end;
    % run 3dDeconvolve
    [parms,errcode] = run_3dDeconv(parms,i);
    if errcode, return; end;
    if parms.vox_flag
      % resample output to T1 resolution
      if parms.resamp_flag
        fprintf('%s: resampling output to T1 resolution...\n',mfilename);
        parms = resamp_output(parms,i);
      end;
      % calculate averages for subcortical ROIs
      if parms.aseg_flag
        fprintf('%s: aseg ROI analysis...\n',mfilename);
        for sem_flag=0:1
          aseg_analysis(parms,i,sem_flag);
        end;
      end;
      % paint voxelwise results to surface
      if parms.paint_flag
        fprintf('%s: resampling output to cortical surface...\n',mfilename);
        parms = paint_output(parms,i);
        % calculate averages for surface ROIs
        if parms.aparc_flag 
          fprintf('%s: aparc ROI analysis...\n',mfilename);
          for sem_flag=0:1
            for j=1:length(parms.aparc_infix)
              aparc_analysis(parms,i,sem_flag,parms.aparc_infix{j});
            end;
          end;
        end;
        if parms.fparc_flag
          fprintf('%s: fparc ROI analysis...\n',mfilename);
          for sem_flag=0:1
            aparc_analysis(parms,i,sem_flag,'fparc');
          end;
        end;
        if parms.ico_presmooth
          fprintf('%s: smoothing output...\n',mfilename);
          parms = smooth_output(parms);
        end;
        if parms.ico_flag
          fprintf('%s: resampling output to icosahedral sphere with order = %d...\n',...
            mfilename,parms.ico_order);
          sphere_output(parms);
        end;
        % create output files for each condition
        fprintf('%s: splitting output...\n',mfilename);
        split_output(parms,i);
      end;
      % threhold coefs with tstats
      if parms.pthresh~=0 && parms.pthresh~=1
        fprintf('%s: thresholding output...\n',mfilename);
        threshold_output(parms);
      end;
      % render surface map images
      if parms.render_flag && parms.paint_flag
        fprintf('%s: rendering output...\n',mfilename);
        render_output(parms,i);
      end;
      % adding voxel-wise outputs 
      if parms.mni_flag
        fprintf('%s: writing mni outputs...\n',mfilename);
        mni_outputs(parms);
      end;
    end;
    % collect results for ROI-average analyses
    if parms.roi_flag
      if parms.aseg_flag
        fprintf('%s: compiling aseg ROI analysis...\n',mfilename);
        for sem_flag=0:1
         roi_analysis(parms,i,sem_flag,'aseg'); 
        end;
      end;
      if parms.aparc_flag
        fprintf('%s: compiling aparc ROI analysis...\n',mfilename);
        for sem_flag=0:1
          for j=1:length(parms.aparc_infix)
            roi_analysis(parms,i,sem_flag,parms.aparc_infix{j});
          end;
        end;
      end;
      if parms.fparc_flag
        fprintf('%s: compiling fparc ROI analysis...\n',mfilename);
        for sem_flag=0:1
          roi_analysis(parms,i,sem_flag,'fparc'); 
        end;
      end;
    end;
    % save scan info
    fname_info = sprintf('%s/%s_info.mat',parms.outdirs{i},parms.outstems{i});
    scan_info = [];
    if ~exist(fname_info,'file') || parms.forceflag
      fprintf('%s: saving scan info to %s...\n',mfilename,fname_info);
      scan_info = get_scan_info(parms,s);
      if ~isempty(scan_info)
        save(fname_info,'scan_info');
      end;
    end;
    % save motion stats
    if parms.mc_flag
      fname_motion = sprintf('%s/%s_motion.mat',...
        parms.outdirs{i},parms.outstems{i});
      if ~exist(fname_motion,'file') || parms.forceflag
        fprintf('%s: saving motion stats to %s...\n',mfilename,fname_motion);
        if isempty(scan_info) && exist(fname_info,'file')
          load(fname_info);
          numTRs = scan_info.numTRs;
        else
          [~,ind_snum] = intersect(parms.snums,s);
          numTRs = parms.nreps(ind_snum);
        end;
        motion_stats = calc_motion_stats(parms,numTRs);
        if ~isempty(motion_stats)
          save(fname_motion,'-struct','motion_stats');
        end;
      end;
    end;
    % save dof's
    fname_dof = sprintf('%s/%s_dof.mat',parms.outdirs{i},parms.outstems{i});
    if ~exist(fname_dof,'file') || parms.forceflag
      fprintf('%s: saving degrees of freedom to %s...\n',mfilename,fname_dof);
      stats_info = load_stats_info(parms);
      if ~isempty(stats_info)
        ind = select_coefs(parms);
        nconds = length(ind);
        dof1 = zeros(nconds,1);
        dof2 = zeros(nconds,1);
        for j=1:nconds
          k = ind(j);
          dof(j) = stats_info(k+1).dofs;
        end;
        condnames = parms.condnames;
        save(fname_dof,'nconds','condnames','dof');
      end;
    end;
    % censor outliers 
    if parms.outliers_flag
      fprintf('%s: censoring outliers...\n',mfilename);  
      censor_outliers(parms,i); 
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function censor_outliers(parms,run) 
  [input,output] = get_outliers_files(parms);
  fname_out = sprintf('%s/%s_censored%0.1f.mat',...
    parms.outdirs{run},parms.outstems{run},parms.outliers_thresh); 
  if ~exist(fname_out,'file')  
    sem_data = []; do_censoring = 0; 
    for i=1:length(input)
      [vol,~,~,volsz] = fs_load_mgh(input{i});
      vsi = squeeze(vol); 
      sem_ind = select_coefs(parms,input{i},volsz) + 1;
      sem_data = [sem_data; vsi(:,sem_ind)];
    end;
    sem_rms = sqrt(mean(sem_data.^2,1));
    outliers_thresh = parms.outliers_thresh;
    excluded_ind = find(sem_rms>outliers_thresh);
    if ~isempty(excluded_ind),
      do_censoring = 1;
    end;
    condnames = parms.condnames;
    excluded_contrast = {parms.condnames{excluded_ind}};
    censoring_files = input;
    save(fname_out,'outliers_thresh','sem_rms','condnames',...
      'excluded_ind','excluded_contrast','do_censoring',...
      'censoring_files');
  else 
    load(fname_out);
  end;
  % censor ico data vertex-native 
  parms = get_ico(parms);
  for h=1:length(parms.hemilist)
    ico_input_file = sprintf('%s/%s-%s_3dDeconv-%s%s',...
      parms.outdir,parms.outstem,parms.ico_stem,parms.hemilist{h},parms.out_ext);
    ico_output_file = sprintf('%s/%s-%s_3dDeconv_censored%0.1f-%s%s',...
      parms.outdir,parms.outstem,parms.ico_stem,parms.outliers_thresh,...
      parms.hemilist{h},parms.out_ext);
    ico_input_file_loc = sprintf('%s-%s_3dDeconv-%s%s',...
      parms.outstem,parms.ico_stem,parms.hemilist{h},parms.out_ext);
    ico_output_file_loc = sprintf('%s-%s_3dDeconv_censored%0.1f-%s%s',...
      parms.outstem,parms.ico_stem,parms.outliers_thresh,...
      parms.hemilist{h},parms.out_ext);
    if ~exist(ico_output_file)
      if exist(ico_input_file)
        if do_censoring
          [vol,M,mr_parms,volsz] = fs_load_mgh(ico_input_file);
          vol_censored = vol;
          % NOTE: vol has interleaved beta and SEM
          ind_excl = [2*excluded_ind,1+2*excluded_ind];
          vol_censored(:,1,1,ind_excl) = nan;
          fs_save_mgh(vol_censored,ico_output_file,M,mr_parms);
        else
          cmd = sprintf('cd %s ;ln -s %s %s ; cd -',parms.outdir,...
            ico_input_file_loc,ico_output_file_loc);
          status = unix(cmd);
          if status,
            error('censor_outliers: failed link to %s\n',ico_input_file);
          end;
        end;
      end;
    end;
  end;
  % censor all roi-based available
  for i=1:length(output)
    [fpath,fstem,fext] = fileparts(output{i});
    fname_out = sprintf('%s/%s_censored%0.1f%s',fpath,fstem,...
      parms.outliers_thresh,fext);
    fname_out_loc = sprintf('%s_censored%0.1f%s',fstem,...
      parms.outliers_thresh,fext);
    if ~exist(fname_out,'file') 
      if do_censoring
        roi_data = []; load(output{i});
        if ~exist('roi_data','var')
          error('censor_outliers: no roi_data for %s\n',output{i});
        end; 
        fields = fieldnames(roi_data);
        for f=1:length(fields)
          for r=1:length(roi_data)
            if length(roi_data(r).(fields{f}))==length(condnames) && ...
                  ~ischar(roi_data(r).(fields{f}))
              roi_data(r).(fields{f})(excluded_ind) = nan;
            end;
          end;
        end;  
        save(fname_out,'roi_data','-v7.3');
      else
        cmd = sprintf('cd %s ; ln -s %s%s %s ; cd -',fpath,fstem,fext,fname_out_loc);
        status = unix(cmd);
        if status,
          error('censor_outliers: failed link to %s\n',output{i});
        end;
      end;
    end;
  end;
  % censor mni
  if parms.mni_flag
    conds = parms.condnames; 
    for i=1:length(conds)
      for t=1:length(parms.stats)
        fname_mni_in = sprintf('%s/contrasts_mni/%s_%s_%s_reg2mni.nii.gz',...
          parms.outdirs{run},parms.taskname,conds{i},parms.stats{t});
        fname_mni_in_loc = sprintf('%s_%s_%s_reg2mni.nii.gz',...
          parms.taskname,conds{i},parms.stats{t});
        fname_mni_out = sprintf('%s/contrasts_mni/%s_%s_%s_reg2mni_censored%0.1f.nii.gz',...
          parms.outdirs{run},parms.taskname,conds{i},parms.stats{t},...
          parms.outliers_thresh);
        fname_mni_out_loc = sprintf('%s_%s_%s_reg2mni_censored%0.1f.nii.gz',...
          parms.taskname,conds{i},parms.stats{t},parms.outliers_thresh);
        fname_empty = sprintf('%s/contrasts_mni/%s_empty_reg2mni.nii.gz',...
          parms.outdirs{run},parms.taskname);
        fname_empty_loc = sprintf('%s_empty_reg2mni.nii.gz',parms.taskname);
        if ~exist(fname_mni_out,'file')
          if exist(fname_mni_in,'file')
            if find(strcmp(conds{i},excluded_contrast))
              if ~exist(fname_empty,'file')
                tmp_dir = sprintf('%s/contrasts_mni/',parms.outdirs{run});
                tmp_file = sprintf('%s/contrasts_mni/tmp_%s_empty_reg2mni.mgh',...
                  parms.outdirs{run},parms.taskname);
                fs_mri_convert(fname_mni_in,tmp_file);
                [vol_tmp,M,mr_parms] = fs_load_mgh(tmp_file);
                vol_tmp_nan = vol_tmp;
                vol_tmp_nan(find(abs(vol_tmp)>0)) = nan;
                fs_save_nifti(vol_tmp_nan,fname_empty,M,mr_parms);
                cmd = sprintf('cd %s; 3drefit -space MNI %s; cd -',...
                  tmp_dir,fname_empty);
                [status,result] = unix(cmd);
                if status
                  error('censor_outliers: failed 3drefit on %s\n',fname_empty);
                end;   
                delete(tmp_file);
              end;
              cmd = sprintf('cd %s/contrasts_mni ; ln -s %s %s ; cd -',...
                parms.outdirs{run},fname_empty_loc,fname_mni_out_loc);
              status = unix(cmd);
              if status,
                error('censor_outliers: failed link %s\n',fname_mni_out);
              end;
            else
              if ~exist(fname_mni_out,'file') 
                cmd = sprintf('cd %s/contrasts_mni ; ln -s %s %s ; cd -',...
                  parms.outdirs{run},fname_mni_in_loc,fname_mni_out_loc);
                status = unix(cmd);
                if status,
                  error('censor_outliers: failed link to %s\n',fname_mni_in);
                end;
              end;
            end;
          else
            error('censor_outliers: %s not found \n',fname_mni_in);    
          end;
        end;
      end;
    end;
  end;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [input_array,output_array] = get_outliers_files(parms)
  input_array = {}; output_array = {};
  % censoring based on surface-native vertex data 
  parms = get_ico(parms);
  for h=1:length(parms.hemilist)
    input_file = sprintf('%s/%s-%s_3dDeconv-%s%s',...
      parms.outdir,parms.outstem,parms.ico_stem,parms.hemilist{h},parms.out_ext);
    if exist(input_file,'file')
       input_array{h} = input_file;
    else
      error('censor_outliers: missing input file %s\n',input_file);
    end;
  end;
  if parms.vox_flag
    if parms.aseg_flag
      for s = 1:length(parms.stats)
        output_file = sprintf('%s/%s_3dDeconv_%s_aseg_roi_data.mat',...
          parms.outdir,parms.outstem,parms.stats{s});
        if exist(output_file,'file')
          output_array{end+1} = output_file;
        else
          error('censor_outliers: missing input file %s\n',output_file);
        end;
      end;
    end;
    if parms.aparc_flag
      for s=1:length(parms.stats)    
        for i=1:length(parms.aparc_infix)
          output_file = sprintf('%s/%s_3dDeconv_%s_%s_roi_data.mat',...
            parms.outdir,parms.outstem,parms.stats{s},parms.aparc_infix{i});
          if exist(output_file,'file')
            output_array{end+1} = output_file;
          else
            error('censor_outliers: missing input file %s\n',output_file);
          end;
        end;
      end;
    end;
    if parms.fparc_flag
      for s=1:length(parms.stats)    
        output_file = sprintf('%s/%s_3dDeconv_%s_fparc_roi_data.mat',...
          parms.outdir,parms.outstem,parms.stats{s});
        if exist(output_file,'file')
          output_array{end+1} = output_file;
        else
          error('censor_outliers: missing input file %s\n',output_file);
        end;
      end;
    end;
  end;
  if parms.roi_flag
    if parms.aseg_flag
      for s=1:length(parms.stats)
        output_file = sprintf('%s/%s_ROI_3dDeconv_%s_aseg_roi_data.mat',...
          parms.outdir,parms.outstem,parms.stats{s});
        if exist(output_file,'file')
          output_array{end+1} = output_file;
        else
          error('censor_outliers: missing input file %s\n',output_file);
        end;
      end;
    end;
    if parms.aparc_flag
      for s=1:length(parms.stats)    
        for i=1:length(parms.aparc_infix)
          output_file = sprintf('%s/%s_ROI_3dDeconv_%s_%s_roi_data.mat',...
            parms.outdir,parms.outstem,parms.stats{s},parms.aparc_infix{i});
          if exist(output_file,'file')
            output_array{end+1} = output_file;
          else
            error('censor_outliers: missing input file %s\n',output_file);
          end;
        end;
      end;
    end;
    if parms.fparc_flag
      for s=1:length(parms.stats)
        output_file = sprintf('%s/%s_ROI_3dDeconv_%s_fparc_roi_data.mat',...
          parms.outdir,parms.outstem,parms.stats{s});
        if exist(output_file,'file')
          output_array{end+1} = output_file;
        else
          error('censor_outliers: missing input file %s\n',output_file);
        end;
      end;
    end;
  end;
  if isempty(input_array), 
    error('censor_outliers: no surf data found to compute outliers');
  end;
  if isempty(output_array), 
    error('censor_outliers: no roi data found to censor');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = get_ico(parms)
  parms.ico_stem = []; 
  if parms.ico_presmooth
    parms.ico_stem = sprintf('sm%d-',parms.ico_presmooth);
  end;
  if parms.ico_flag    
    if parms.ico_trunc_flag && parms.ico_order~=7
      parms.ico_stem = sprintf('%sico7-',parms.ico_stem);
    else
      parms.ico_stem = sprintf('%sico%d',parms.ico_stem,parms.ico_order);
    end;
    if parms.ico_postsmooth
      parms.ico_stem = sprintf('%s-sm%d',parms.ico_stem,parms.ico_postsmooth); 
    end;
    if parms.ico_trunc_flag && parms.ico_order~=7
      parms.ico_stem = sprintf('%strunc%d',parms.ico_stem,parms.ico_order);
    end;
  end;
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode] = run_3dDeconv(parms,i)
  errcode = 0;
  % select input files
  s = parms.snum_sets{i};
  [~,d_ind] = intersect(parms.snums,s);
  if length(s)==1
    fname_in = parms.data_fnames{d_ind};
    if parms.mc_flag
      parms.fname_motion = parms.motion_fnames{d_ind};
    else
      parms.fname_motion = [];
    end;
    stim_fnames = parms.stim_fnames{d_ind};
  else
    fname_in = parms.data_fnames(d_ind);
    if parms.mc_flag
      parms.fname_motion = parms.motion_fnames(d_ind);
    else
      parms.fname_motion = [];
    end;
    stim_fnames = parms.stim_fnames(d_ind);
  end;
  % set output
  parms.outdir = parms.outdirs{i};
  parms.outstem = parms.outstems{i};
  parms.regfile = parms.reg_fnames{i};
  % run mmil_3dDeconv
  args = mmil_parms2args(parms,parms.deconv_tags);
  [output,errcode] = mmil_3dDeconv(fname_in,stim_fnames,args{:});
  if errcode 
    fprintf('%s: ERROR: mmil_3dDeconv failed\n',mfilename);
    return;
  end;
  all_fname_out = struct2files(output);
  for i=1:length(all_fname_out)
    if ~exist(all_fname_out{i},'file')
      [~,name,ext] = fileparts(all_fname_out{i});
      error('3dDeconv output file %s%s not found',name,ext);
    end;
  end;
  parms.output = output;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resamp_output(parms,i)
  % set parameters
  s = parms.snum_sets{i};
  [~,s_ind] = intersect(parms.snums,s);
  tp = [];
  tp.forceflag = parms.overwrite_flag;
  tp.fname_regdat = parms.reg_fnames{s_ind(1)};
  tp.fname_ref = sprintf('%s/mri/T1.mgz',parms.fspath);
  if ~exist(tp.fname_ref,'file')
    error('file %s not found',tp.fname_ref);
  end;
  % first frame of BOLD data
  fname_reg = parms.data_fnames{s_ind(1)};
  [fpath,fstem,fext] = fileparts(fname_reg);
  if ~isempty(parms.out_ext)
    fext = parms.out_ext;
  end;
  tp.fname_out = sprintf('%s/%s_resT1_f0%s',fpath,fstem,fext);
  tp.frames = 1;
  args = mmil_parms2args(tp);
  [M_ref2reg,subj,inplane,slicethick] = ...
    mmil_resample_by_regdat(fname_reg,args{:});
  fname_reg = parms.output.vox_glm_results;
  [fpath,fstem,fext] = fileparts(fname_reg);
  if ~isempty(parms.out_ext)
    fext = parms.out_ext;
  end;
  tp.fname_out = sprintf('%s/%s_resT1%s',fpath,fstem,fext);
  tp.frames = [];
  args = mmil_parms2args(tp);
  [M_ref2reg,subj,inplane,slicethick] = ...
    mmil_resample_by_regdat(fname_reg,args{:});
  parms.output.vox_glm_results = tp.fname_out;
  % iresp
  if parms.iresp_flag
    for c=1:length(parms.output.iresp_mgz)
      fname_reg = parms.output.iresp_mgz{c};
      [fpath,fstem,fext] = fileparts(fname_reg);
      tp.fname_out = sprintf('%s/%s_resT1%s',fpath,fstem,fext);
      tp.frames = [];
      args = mmil_parms2args(tp);
      [M_ref2reg,subj,inplane,slicethick] = ...
        mmil_resample_by_regdat(fname_reg,args{:});
      parms.output.iresp_mgz{c} = tp.fname_out;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = paint_output(parms,i)
  s = parms.snum_sets{i};
  [~,s_ind] = intersect(parms.snums,s);
  if parms.resamp_flag
    fname_regdat = [];
  else
    fname_regdat = parms.reg_fnames{s_ind(1)};
  end;
  parms.outtype = regexprep(parms.out_ext,'\.','');
  parms.regfile = fname_regdat;
  args = mmil_parms2args(parms,parms.paint_tags);
  % paint GLM stats
  fname_in = parms.output.vox_glm_results;
  parms.output.vox_glm_results = fs_paint(parms.subj,fname_in,args{:});
  % paint impulse response functions
  if parms.iresp_flag
    [fpath,fstem,fext] = fileparts(fname_in);
    if ~isempty(parms.out_ext)
      fext = parms.out_ext;
    end;
    for c=1:length(parms.output.iresp_mgz)
      fname_in = parms.output.iresp_mgz{c};
      parms.output.iresp_mgz{c} = fs_paint(parms.subj,fname_in,args{:});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = smooth_output(parms)
  for f=1:length(parms.output.vox_glm_results)
    fname_in = parms.output.vox_glm_results{f};
    [fpath,fstem,fext] = fileparts(fname_in);
    n = regexp(fstem,'(?<stem>.+)-(?<hemi>\wh)','names');
    if ~isempty(n)
      hemi = n.hemi;
      fstem = n.stem;
    else
      hemi = [];
    end;
    fname_out = sprintf('%s/%s-sm%d-%s%s',...
      fpath,fstem,parms.ico_presmooth,hemi,fext);
    parms.output.vox_glm_results{f} = fs_surf2surf(fname_in,parms.subj,...
      'fname_out',fname_out,...
      'smooth_out',parms.ico_presmooth,...
      'forceflag',parms.forceflag);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sphere_output(parms)
  for f=1:length(parms.output.vox_glm_results)
    fname_in = parms.output.vox_glm_results{f};
    [fpath,fstem,fext] = fileparts(fname_in);
    n = regexp(fstem,'(?<stem>.+)-(?<hemi>\wh)','names');
    if ~isempty(n)
      hemi = n.hemi;
      fstem = n.stem;
    else
      hemi = [];
    end;
    if parms.ico_trunc_flag && parms.ico_order~=7
      fstem = sprintf('%s-ico7',fstem);
    else
      fstem = sprintf('%s-ico%d',fstem,parms.ico_order);
    end;
    if parms.ico_postsmooth
      fstem = sprintf('%s-sm%d',fstem,parms.ico_postsmooth);
    end;
    fname_out = sprintf('%s/%s-%s%s',fpath,fstem,hemi,fext);
    tparms = [];
    tparms.fname_out = fname_out;
    tparms.trgsubj = 'ico';
    tparms.hemi = hemi;
    if ~parms.ico_trunc_flag
      tparms.icolevel = parms.ico_order;
    else
      tparms.icolevel = 7;
    end;
    tparms.smooth_out = parms.ico_postsmooth;
    tparms.subjdir = parms.subjdir;
    tparms.intype = parms.outtype;
    tparms.forceflag = parms.forceflag;
    tparms.verbose = 0;
    args = mmil_parms2args(tparms);
    fname_out = fs_surf2surf(fname_in,parms.subj,args{:});
    % create truncated copy with ico_order number of vertices
    if parms.ico_trunc_flag && parms.ico_order~=7
      fname_in = fname_out;
      fstem = sprintf('%s-trunc%d',fstem,parms.ico_order);
      fname_out = sprintf('%s/%s-%s%s',fpath,fstem,hemi,fext);
      if ~exist(fname_out,'file') || parms.forceflag
        nverts = parms.ico_nverts(parms.ico_order);
        vals = fs_load_mgh(fname_in);
        vals = vals(1:nverts,:,:,:);
        fs_save_mgh(vals,fname_out);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stats_info = load_stats_info(parms)
  stats_info = [];
  if ~isfield(parms.output,'info')
    error('no stats info');
  end;
  if exist(parms.output.info,'file')
    load(parms.output.info);
    if isempty(stats_info)
      error('file %s did not contain expected stats info',parms.output.info);
    end;
  else
    sprintf('WARNING: stats_info file not found\n');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshold_output(parms)
  stats_info = load_stats_info(parms);
  parms.fname_stats = parms.output.vox_glm_results;
  if ~iscell(parms.fname_stats), 
    parms.fname_stats = {parms.fname_stats}; 
  end;
  for f=1:length(parms.fname_stats)
    fname_in = parms.fname_stats{f};
    [fpath,fstem,fext] = fileparts(fname_in);
    n = regexp(fstem,'(?<stem>.+)-(?<hemi>\wh)','names');
    if ~isempty(n)
      hemi = n.hemi;
      fstem = n.stem;
    else
      hemi = [];
    end;
    if ~isempty(hemi)
      fname_in = sprintf('%s/%s-%s%s',...
        fpath,fstem,hemi,fext);
      fname_out = sprintf('%s/%s-pthresh%0.1e-%s%s',...
        fpath,fstem,parms.pthresh,hemi,fext);
    else
      fname_in = sprintf('%s/%s%s',...
        fpath,fstem,fext);
      fname_out = sprintf('%s/%s-pthresh%0.1e%s',...
        fpath,fstem,parms.pthresh,fext);
    end;
    if ~exist(fname_out,'file') || parms.forceflag
      [vol,M] = fs_load_mgh(fname_in,[],[],0,1);
      nframes = size(vol,4);
      % select coef frames (all even)
      ind = select_coefs(parms,fname_in,nframes);
      vol_out = zeros([size(vol,1),size(vol,2),size(vol,3),length(ind)]);
      for j=1:length(ind)
        k = ind(j);
        if k+1>nframes
          error('expected more frames in %s',fname_in);
        end;
        tmp_vol_coef = vol(:,:,:,k);
        % can make assumptions about where tstats are relative to coefs
        tmp_vol_sem =  vol(:,:,:,k+1);
        tmp_vol_tstat = tmp_vol_coef./tmp_vol_sem;
        dofs = stats_info(k+1).dofs;
        tthresh = abs(tinv(1-parms.pthresh,dofs));
        tmp_vol_coef(abs(tmp_vol_tstat)<tthresh)=0;
        vol_out(:,:,:,j) = tmp_vol_coef;
      end;
      % save thresholded coefs in multi-frame file
      fs_save_mgh(vol_out,fname_out,M);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_coef = select_coefs(parms,fname_in,nframes)
  if ~exist('fname_in','var'), fname_in = []; end;
  if ~exist('nframes','var'), nframes = []; end;
  % check that number of frames matches twice the number of contrasts plus 1
  if ~isempty(fname_in) && ~isempty(nframes)
    nexpect = 1 + 2*parms.nconds;
    if nframes ~= nexpect
      error('mismatch between number of frames in %s (%d) and expected (%d)',...
        fname_in,nframes,nexpect);
    end;
  end;
  % select coef frames (all even)
  ind_coef = 1 + [1:2:2*parms.nconds];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [contdir,fstem] = split_output(parms,i,subj_flag,post_thresh_flag)
  if ~exist('subj_flag','var') || isempty(subj_flag), subj_flag = 0; end;
  if ~exist('post_thresh_flag','var') || isempty(post_thresh_flag), post_thresh_flag = 0; end;
  indir = parms.outdirs{i};
  contdir = [parms.outdir '/contrasts'];
  fstem = set_fstem(parms,i,subj_flag,post_thresh_flag);
  mmil_mkdir(contdir);
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    fname_in = sprintf('%s/%s-%s%s',indir,fstem,hemi,parms.out_ext);
    if ~exist(fname_in,'file')
      error('file %s not found',fname_in);
    end;
    vol = [];
    if ~post_thresh_flag || ismember(parms.pthresh,[0,1])
      for t = 1:length(parms.stats)
        stat = parms.stats{t};    
        for k=1:parms.nconds
          condname = parms.condnames{k};
          fname = sprintf('%s/%s_%s_%s-%s.mgz',contdir,fstem,condname,stat,hemi);
          if ~exist(fname,'file') || parms.forceflag
            if isempty(vol)
              vol = fs_load_mgh(fname_in);
            end
            ind = select_coefs(parms,fname_in,size(vol,4)) + (t-1);
            vol_stat = vol(:,:,:,ind);
            nframes = size(vol_stat,4);
            if nframes~=parms.nconds
              error('mismatch between number of stimuli (%d) plus contrasts (%d) and number of frames (%d) in %s',...
                parms.nstims,parms.ncontrasts,nframes,fname_in);
            end;
            fs_save_mgh(vol_stat(:,:,:,k),fname);
          end;
        end;
      end;    
    else
      for k=1:parms.nconds
        condname = parms.condnames{k};
        fname = sprintf('%s/%s_%s-%s.mgz',contdir,fstem,condname,hemi);
        if ~exist(fname,'file') || parms.forceflag
          if isempty(vol)
            vol = fs_load_mgh(fname_in);
            nframes = size(vol,4);
            if nframes~=parms.nconds
              error('mismatch between number of stimuli (%d) plus contrasts (%d) and number of frames (%d) in %s',...
                parms.nstims,parms.ncontrasts,nframes,fname_in);
            end;
          end;
          fs_save_mgh(vol(:,:,:,k),fname);
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = set_fstem(parms,i,subj_flag,post_thresh_flag)
  fstem = [parms.outstems{i} '_3dDeconv'];
  if parms.resamp_flag
    fstem = [fstem '_resT1'];
  end;
  if parms.ico_presmooth
    fstem = sprintf('%s-sm%d',fstem,parms.ico_presmooth);
  end;
  if ~subj_flag && parms.ico_flag
    if parms.ico_trunc_flag && parms.ico_order~=7
      fstem = sprintf('%s-ico7',fstem);
    else
      fstem = sprintf('%s-ico%d',fstem,parms.ico_order);
    end;
    if parms.ico_postsmooth
      fstem = sprintf('%s-sm%d',fstem,parms.ico_postsmooth);
    end;
    if parms.ico_trunc_flag && parms.ico_order~=7
      fstem = sprintf('%s-trunc%d',fstem,parms.ico_order);
    end;
  end;
  if post_thresh_flag && ~ismember(parms.pthresh,[0,1])
    fstem = sprintf('%s-pthresh%0.1e',fstem,parms.pthresh);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aseg_analysis(parms,i,sem_flag)
  fname_in = parms.output.vox_glm_results;  
  [M,volsz] = mmil_load_mgh_info(fname_in,parms.forceflag,parms.outdirs{i});
  nframes = volsz(4);
  parms.outdir = parms.outdirs{i};
  parms.outstem = [parms.outstems{i} '_3dDeconv'];
  if sem_flag
    parms.outstem = [parms.outstem '_sem'];
  else
    parms.outstem = [parms.outstem '_beta'];
  end;
  parms.frames = select_coefs(parms,fname_in,nframes);
  if sem_flag
    parms.frames = parms.frames + 1;
  end;
  if parms.resamp_flag
    parms.M_reg = [];
  else
    parms.M_reg = load_regdat(parms,i);
  end;
  parms.verbose = (parms.verbose>0);
  parms.output = [];
  args = mmil_parms2args(parms,parms.aseg_tags);
  mmil_aseg_analysis(fname_in,parms.fspath,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = roi_analysis(parms,i,sem_flag,seg_flag)
  errcode = 0;
  if strcmp(seg_flag,'aseg')
    fname_in = parms.output.roi_vol_glm_results;
    if ~iscell(fname_in), fname_in = {fname_in}; end;
  elseif strcmp(seg_flag,'aparc') || strcmp(seg_flag,'aparc.a2009s')
    ind = find(strcmp(parms.aparc_names,seg_flag));
    fname_in = parms.output.roi_surf_glm_results(ind);
  elseif strcmp(seg_flag,'fparc')
    num_fparc = length(parms.fnames_fparc_resamp);
    fname_in = cell(num_fparc,1);
    for f=1:num_fparc
      ind = find(strcmp(parms.fnames_aparc,parms.fnames_fparc_resamp{f}));
      fname_in(f) = parms.output.roi_surf_glm_results(ind);
    end;
  else
    errcode = 1; return;
  end;
  parms.outdir = parms.outdirs{i};
  parms.outstem = [parms.outstems{i} '_ROI_3dDeconv'];
  if sem_flag
    parms.outstem = [parms.outstem '_sem'];
  else
    parms.outstem = [parms.outstem '_beta'];
  end;
  fname_out = sprintf('%s/%s_%s_roi_data.mat',...
    parms.outdir,parms.outstem,seg_flag);
  if ~exist(fname_out,'file')
    roi_data = struct; pos = 0;
    for h=1:length(fname_in)
      if exist(fname_in{h},'file')
        required_tags = {'roinames','roicodes',...
          'nvals','nvals_invalid','nvals_valid'};
        tmp_info = mmil_load_matfile(fname_in{h},required_tags);
      else
        errcode = 1; return;
      end;
      nframes = size(tmp_info.vol,2); 
      parms.frames = select_coefs(parms,fname_in,nframes);
      if sem_flag
        parms.frames = parms.frames + 1;
      end;
      for j=1:length(tmp_info.roinames)
        i = pos+j;
        roi_data(i).roicode = tmp_info.roicodes(j);
        roi_data(i).roiname = tmp_info.roinames{j};
        roi_data(i).vals = [];
        roi_data(i).weights = [];
        roi_data(i).nvals = tmp_info.nvals(j);
        roi_data(i).nvals_valid = tmp_info.nvals_valid(j);
        roi_data(i).nvals_invalid = tmp_info.nvals_invalid(j);
        roi_data(i).avg = tmp_info.vol(j,parms.frames);
        roi_data(i).median = [];
        roi_data(i).stdv = [];
      end;
      pos = length(tmp_info.roinames);
    end;
    save(fname_out,'roi_data','-v7.3');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M_reg = load_regdat(parms,i)
  % get scan index
  s = parms.snum_sets{i};
  [tmp,s_ind] = intersect(parms.snums,s);
  s_ind = s_ind(1);
  % get info for T1 reference image
  fname_ref = sprintf('%s/mri/T1.mgz',parms.fspath);
  if ~exist(fname_ref,'file')
    error('file %s not found',fname_ref);
  end;
  [M_ref,nvox_ref] = mmil_load_mgh_info(fname_ref,parms.forceflag);
  % get info for BOLD data
  fname_reg = parms.data_fnames{s_ind};
  [M_reg,nvox_reg] = mmil_load_mgh_info(fname_reg,parms.forceflag);
  % load regdat file
  fname_regdat = parms.reg_fnames{s_ind};
  M_reg = fs_read_regdat(fname_regdat,...
    'tk2ras_flag',1,...
    'M_ref',M_ref,...
    'M_reg',M_reg,...
    'nvox_ref',nvox_ref,...
    'nvox_reg',nvox_reg);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aparc_analysis(parms,i,sem_flag,annot_name)
  if ~exist('annot_name','var') || isempty(annot_name)
    annot_name = 'aparc';
  end;
  outstem = sprintf('%s/%s_3dDeconv',parms.outdirs{i},parms.outstems{i});
  if sem_flag
    outstem = [outstem '_sem'];
  else
    outstem = [outstem '_beta'];
  end;
  fname_out = [outstem '_' annot_name '_roi_data.mat'];
  if ~exist(fname_out,'file') || parms.forceflag
    roi_data = [];
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      fname_in =  parms.output.vox_glm_results{h};
      [~,volsz] = mmil_load_mgh_info(fname_in,parms.forceflag,parms.outdirs{i});
      nframes = volsz(4);
      parms.hemi = hemi;
      parms.annot_name = annot_name;
      parms.fname_aparc = set_fname_aparc(parms,hemi,annot_name);
      parms.frames = select_coefs(parms,fname_in,nframes);
      parms.verbose = (parms.verbose > 0);
      if sem_flag
        parms.frames = parms.frames + 1;
      end;
      args = mmil_parms2args(parms,parms.surf_roi_tags);
      tmp_roi_data = mmil_surf_roi(fname_in,args{:});
      roi_data = [roi_data,tmp_roi_data];
    end;
    save(fname_out,'roi_data');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function render_output(parms,i)
 [contdir,fstem] = split_output(parms,i,1,1);
  parms.surfdir = [parms.outdir '/surfs'];
  parms.outdir = [parms.outdir '/plots'];
  parms.surfname = parms.view_surfname;
  mmil_mkdir(parms.outdir);
  for h=1:length(parms.hemilist)
    hemi = parms.hemilist{h};
    for k=1:parms.nconds
      condname = parms.condnames{k};
      fname = sprintf('%s/%s_%s-%s.mgz',contdir,fstem,condname,hemi);
      [fpath,instem,fext] = fileparts(fname);
      for v=1:length(parms.viewlist)
        parms.view = parms.viewlist{v};
        parms.outstem = sprintf('%s-%s-%s',...
          instem,parms.surfname,parms.view);
        fname_out = [parms.outdir '/' parms.outstem '.tif'];
        if ~exist(fname_out,'file') || parms.forceflag
          args = mmil_parms2args(parms,parms.sview_tags);
          sv_surf_view(parms.subj,fname,args{:});
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion_stats = calc_motion_stats(parms,numTRs)
  motion_stats = [];
  motion_data = zeros(0,6);
  max_dz = 0; max_dy = 0; max_dx = 0;
  max_rz = 0; max_rx = 0; max_ry = 0;
  mean_trans = 0; max_trans = 0; mean_rot = 0; max_rot = 0;
  mean_motion = 0; mode_motion = 0; med_motion = 0;
  min_motion = Inf; max_motion = 0; mean_accel = 0;
  subthresh_nvols = 0; subthresh_perc = 0;
  thresholds = 0; nvols = 0;
  if ~iscell(parms.fname_motion)
    fnames = {parms.fname_motion};
  else
    fnames = parms.fname_motion;
  end;
  nscans = length(fnames);
  for i=1:nscans
    fname_motion = fnames{i};
%    keyboard
    tmp_data_orig = mmil_load_motion_1D(fname_motion,...
      'skipTRs',parms.skipTRs,'nframes',numTRs(i));
    % calculate difference in head position or rotation
    [motion_stats_orig,motion_fd_orig] = mmil_motion_stats(tmp_data_orig,...
      parms.motion_radius,parms.motion_absflag,...
      parms.censor_thresh,parms.motion_nodyflag);
    % resp filter motion estimates to remove respiration noise
    if parms.mc_resp_flag
      tmp_data = abcd_resp_filter_motion(tmp_data_orig,parms.TR,...
        parms.resp_low,parms.resp_high);
      % calculate difference in head position or rotation
      [motion_stats,motion_fd] = mmil_motion_stats(tmp_data,...
        parms.motion_radius,parms.motion_absflag,...
        parms.censor_thresh,parms.motion_nodyflag); %,motion_fd_orig);
    else
      tmp_data = tmp_data_orig;
      motion_stats = motion_stats_orig;
      motion_fd = motion_fd_orig;
    end;
    % combine across scans      
    max_dx = max(max_dx,motion_stats.max_dx);
    max_dy = max(max_dy,motion_stats.max_dy);
    max_dz = max(max_dz,motion_stats.max_dz);
    max_rx = max(max_rx,motion_stats.max_rx);
    max_ry = max(max_ry,motion_stats.max_ry);
    max_rz = max(max_rz,motion_stats.max_rz);
    max_trans = max(max_trans,motion_stats.max_trans);
    max_rot = max(max_rot,motion_stats.max_rot);
    max_motion = max(max_motion,motion_stats.max_motion);
    min_motion = min(min_motion,motion_stats.min_motion);
    mean_trans = mean_trans + motion_stats.mean_trans;
    mean_rot = mean_rot + motion_stats.mean_rot;
    mean_motion = mean_motion + motion_stats.mean_motion;
    mode_motion = mode_motion + motion_stats.mode_motion;
    med_motion = med_motion + motion_stats.med_motion;
    mean_accel = mean_accel + motion_stats.mean_accel;
    subthresh_nvols = subthresh_nvols + motion_stats.subthresh_nvols;
    subthresh_perc = subthresh_perc + motion_stats.subthresh_perc;
    thresholds = motion_stats.thresholds;
    nvols = nvols + motion_stats.nvols;
    motion_data = cat(1,motion_data,tmp_data);
  end;
  if nscans > 1
    mean_trans = mean_trans / nscans;
    mean_rot = mean_rot / nscans;
    mean_motion = mean_motion / nscans;
    mode_motion = mode_motion / nscans;
    med_motion = med_motion / nscans;
    mean_accel = mean_accel / nscans;
    subthresh_perc = subthresh_perc / nscans;
  end;
  motion_stats = struct(...
    'motion_data',motion_data,...
    'max_dx',max_dx,...
    'max_dy',max_dy,...
    'max_dz',max_dz,...
    'max_rx',max_rx,...
    'max_ry',max_ry,...
    'max_rz',max_rz,...
    'mean_trans',mean_trans,...
    'max_trans',max_trans,...
    'mean_rot',mean_rot,...
    'max_rot',max_rot,...
    'mean_motion',mean_motion,...
    'mode_motion',mode_motion,...
    'med_motion',med_motion,...
    'min_motion',min_motion,...
    'max_motion',max_motion,...
    'mean_accel',mean_accel,...
    'subthresh_nvols',subthresh_nvols,...
    'subthresh_perc',subthresh_perc,...
    'thresholds',thresholds,...
    'nvols',nvols);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scan_info = get_scan_info(parms,snums)
  nscans = length(snums);
  scan_info = [];
  scan_info.nscans = nscans;
  scan_info.snums = snums;
  scan_info.ScanInfo = parms.ScanInfo;
  scan_info.SessInfo = parms.SessInfo;
  scan_info.fnames_motion = parms.fname_motion;
  scan_info.fnames_data = {};
  for i=1:nscans
    s = snums(i);
    [~,j] = intersect(parms.snums,s);
    fname_data = parms.data_fnames{j};
    scan_info.fnames_data{i} = fname_data;
    scan_info.TRs(i) = parms.ScanInfo(s).TR/1000;
    scan_info.nreps(i) = parms.ScanInfo(s).nreps;
    % NOTE: nreps and numTRs may be identical
    %       except for APE sequence, where numTRs = nreps/2          
    [M,volsz] = mmil_load_mgh_info(fname_data,parms.forceflag,parms.outdir);
    scan_info.numTRs(i) = volsz(4);
    % adjust TR if numTRs == nreps/2 (i.e. for APE sequence)
    if scan_info.numTRs(i) == scan_info.nreps(i)/2
      scan_info.TRs(i) = scan_info.TRs(i) * scan_info.nreps(i) / scan_info.numTRs(i);
    end;
  end;
  if nscans==1
    scan_info.fnames_data = scan_info.fnames_data{1};
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resamp_fparc(parms)
  if parms.fparc_flag
    if ~isempty(parms.fnames_fparc)
      if parms.verbose>1
        fprintf('%s: resample_fparc_annot...\n',mfilename);
      end;
      parms = resample_fparc_annot(parms);
    else
      parms.fparc_flag = 0;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mni_outputs(parms)
  tp = [];
  tp.tasknames = parms.taskname;  
  tp.infix = parms.infix; 
  tp.outdir = sprintf('%s/contrasts_mni',parms.outdir);
  tp.concat_flag = 0; 
  tp.smooth_fwhm = 5;
  tp.forceflag = parms.forceflag;
  tp.analysis_outfix = parms.analysis_outfix;
  args = mmil_parms2args(tp);
  abcd_reg2mni_taskBOLD(parms.cpath,parms.fspath,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_fname_out = struct2files(output,varargin)
  nomgz_flag = 0;
  for k=1:length(varargin)
    if strcmpi(varargin{k},'nomgz'), nomgz_flag = 1; end;
  end;
  all_fname_out = {};
  output_fields = fieldnames(output);
  for i=1:length(output_fields)
    fname_out_data = output.(output_fields{i});
    if iscell(fname_out_data)
      for j=1:length(fname_out_data)
        all_fname_out{end+1} = fname_out_data{j};    
      end;
    else
     all_fname_out{end+1} = fname_out_data;        
    end;
  end;
  if nomgz_flag == 1  
    all_fname_out_tmp = {};  
    for i=1:length(all_fname_out)
      [~,~,ext] = fileparts(all_fname_out{i});
      if ~strcmp(ext,'.mgz')
        all_fname_out_tmp{end+1} = all_fname_out{i};  
      end;
    end;
    all_fname_out = all_fname_out_tmp;
  end;
  all_fname_out = all_fname_out';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

