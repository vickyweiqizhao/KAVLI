library(R.matlab)
library(readxl)
library(dplyr)
library(tidyr)

####################### paths and importing functions
datapath = '/space/abcd-sync/4.0/tabulated/released/'
funcpath = '/home/wez025/matlab/ABCD_code/github/cmig_library/'
outpath = '/home/wez025/matlab/ABCD_release_3.0/Released/Data/'
source(paste0(funcpath, 'loadtxt.R'))

formattxt <- function(filename, datapath){
    tmp = read.delim(paste0(datapath, filename), skip=2, header=F, sep = "\t")
    cols = read.delim(paste0(datapath, filename), header=T, sep = "\t")
    colnames(tmp) = colnames(cols)
    tmp = tmp[, -c(which(colnames(tmp) %in% c('collection_id','dataset_id')))]
    return(tmp)
}
######################################################################################################
####################### COMBINING: task fmri degrees of freedom & qc flags
######################################################################################################

dat = read.csv('/home/wez025/matlab/ABCD_release_3.0/Released/Data/qcflags_dof_nts_meanmotion_release4.0.csv')
dat = dat[,c("src_subject_id", 'interview_date', 'interview_age', 'eventname', 'mri_info_visitid', 'imgincl_nback_include')]
dat = dat[which(dat$imgincl_nback_include == 1),]
dat = dat[which(dat$eventname == 'baseline_year_1_arm_1'),]

write.csv(dat, paste0('/home/wez025/matlab/ABCD_code/KAVLI/data/', 'subjects_passQC_4.0.csv'), row.names = FALSE)

# fixed_vars = c('src_subject_id','eventname','interview_date', 'interview_age', 'sex')

# imginc = formattxt('abcd_imgincl01.txt', datapath) ## N = 17485
# imginc = imginc[, -which(colnames(imginc) %in% c('visit', 'abcd_imgincl01_id', 'collection_title', 'study_cohort_name'))]

# mri01 = formattxt('abcd_mri01.txt', datapath) ## N = 17485
# mri01 = mri01[,c(fixed_vars, 'mri_info_visitid')]

# fsqc = formattxt('freesqc01.txt', '/home/wez025/matlab/ABCD_release_3.0/released/') ## N = 17438
# fsqc = fsqc[,c(fixed_vars,'fsqc_qc')]
# fsqc$fsqc_qc[fsqc$fsqc_qc==0] = NA

# rsfmri_nts = formattxt('abcd_betnet02.txt', datapath) # N = 16874
# rsfmri_nts = rsfmri_nts[,c(fixed_vars, 'rsfmri_c_ngd_ntpoints')]

# rsfmri_fd = formattxt('mriqcrp10301.txt', datapath) ##N = 17496
# rsfmri_fd = rsfmri_fd[,c(fixed_vars,'iqc_rsfmri_all_mean_motion')]

# nbackdof1 = formattxt('nbackr101.txt', datapath) ## N = 14374
# nbackdof1 = nbackdof1[,c(fixed_vars, 'tfmri_nback_run1_beta_dof')]
# nbackdof2 = formattxt('nbackr201.txt', datapath) ## N = 14073
# nbackdof2 = nbackdof2[,c(fixed_vars, 'tfmri_nback_r2_beta_dof')]

# sstdof1 = formattxt('mrisstr1bw01.txt', datapath) ## N = 14430
# sstdof1 = sstdof1[, c(fixed_vars, 'tfmri_sstr1_beta_dof')]
# sstdof2 = formattxt('mrisstr2bw01.txt', datapath) ## N = 13893
# sstdof2 = sstdof2[,c(fixed_vars, 'tfmri_sstr2_beta_dof')]

# middof1 = formattxt('abcd_midr1bwp102.txt', datapath) ## N = 14650
# middof1 = middof1[,c(fixed_vars, 'tfmri_mid_run1_beta_dof')]
# middof2 = formattxt('midr2bwp102.txt', datapath) ## N = 14290
# middof2 = middof2[,c(fixed_vars, 'tfmri_mr2_beta_dof')]

# tmp = formattxt('mriqcrp10301.txt', datapath) ## N = 17485
# iqc_t1 = tmp[,c(fixed_vars, 'iqc_t1_ok_ser')]
# iqc_t1$iqc_t1_ok_ser[iqc_t1$iqc_t1_ok_ser<=0] = NA

# tmp = formattxt('mriqcrp20301.txt', datapath) ## N = 17485
# tqc_mid = tmp[,c(fixed_vars,'iqc_mid_ok_ser')] ## timeseries passed rawQC
# tqc_mid$iqc_mid_ok_ser[tqc_mid$iqc_mid_ok_ser<=0] = NA

# tmp = formattxt('mriqcrp20301.txt', datapath) ## N = 17485
# tqc_nback = tmp[,c(fixed_vars, 'iqc_nback_ok_ser')]
# tqc_nback$iqc_nback_ok_ser[tqc_nback$iqc_nback_ok_ser<=0] = NA

# tmp = formattxt('mriqcrp20301.txt', datapath) ## N = 17485
# tqc_sst = tmp[,c(fixed_vars,'iqc_sst_ok_ser')]
# tqc_sst$iqc_sst_ok_ser[tqc_sst$iqc_sst_ok_ser<=0] = NA

# tmp = formattxt('abcd_auto_postqc01.txt', datapath) ## N = 17485
# postqc = tmp[,c(fixed_vars, 'apqc_fmri_bounwarp_flag', 'apqc_fmri_regt1_rigid', 'apqc_fmri_fov_cutoff_dorsal', 'apqc_fmri_fov_cutoff_ventral')]
# postqc$apqc_fmri_bounwarp_flag[postqc$apqc_fmri_bounwarp_flag !=1] = NA
# postqc$apqc_fmri_regt1_rigid[postqc$apqc_fmri_regt1_rigid>=19] = NA
# postqc$apqc_fmri_fov_cutoff_dorsal[postqc$apqc_fmri_fov_cutoff_dorsal>=65] = NA
# postqc$apqc_fmri_fov_cutoff_ventral[postqc$apqc_fmri_fov_cutoff_ventral>=60] = NA

# tmp = formattxt('abcd_fmriqc01.txt', datapath) ## N = 17485
# fmri_postqc_qc = tmp[,c(fixed_vars, 'fmri_postqc_qc')]
# fmri_postqc_qc$fmri_postqc_qc[which(is.na(fmri_postqc_qc$fmri_postqc_qc))] = 999
# fmri_postqc_qc$fmri_postqc_qc[fmri_postqc_qc$fmri_postqc_qc==0] = NA

# midqc = Reduce(function(x, y) merge(x, y, all=TRUE), list(iqc_t1, tqc_mid, postqc, fsqc, fmri_postqc_qc))
# midqc$imgincl_taskres_mid_include = 0
# midqc$imgincl_taskres_mid_include[complete.cases(midqc[,grep('iqc_t1_ok_ser', names(midqc)):grep('fmri_postqc_qc', names(midqc))])] = 1
# midqc = midqc[,c(fixed_vars, 'imgincl_taskres_mid_include')]

# sstqc = Reduce(function(x, y) merge(x, y, all=TRUE), list(iqc_t1, tqc_sst, postqc, fsqc, fmri_postqc_qc))
# sstqc$imgincl_taskres_sst_include = 0
# sstqc$imgincl_taskres_sst_include[complete.cases(sstqc[,grep('iqc_t1_ok_ser', names(sstqc)):grep('fmri_postqc_qc', names(sstqc))])] = 1
# sstqc = sstqc[,c(fixed_vars, 'imgincl_taskres_sst_include')]

# nbackqc = Reduce(function(x, y) merge(x, y, all=TRUE), list(iqc_t1, tqc_nback, postqc, fsqc, fmri_postqc_qc))
# nbackqc$imgincl_taskres_nback_include = 0
# nbackqc$imgincl_taskres_nback_include[complete.cases(nbackqc[,grep('iqc_t1_ok_ser', names(nbackqc)):grep('fmri_postqc_qc', names(nbackqc))])] = 1
# nbackqc = nbackqc[,c(fixed_vars, 'imgincl_taskres_nback_include')]

# # qc_flags = Reduce(function(x, y) merge(x, y, by.x = c('src_subject_id', 'eventname', 'interview_date', 'interview_age', 'sex'), by.y = c('src_subject_id', 'eventname', 'interview_date', 'interview_age', 'sex')), list(imginc, mri01,fsqc, rsfmri_nts, rsfmri_fd, nbackdof1, nbackdof2, sstdof1, sstdof2, middof1, middof2))
# qc_flags = Reduce(function(x, y) merge(x, y, all=TRUE), list(imginc, mri01,fsqc, rsfmri_nts, rsfmri_fd, nbackdof1, nbackdof2, sstdof1, sstdof2, middof1, middof2))

# qc_flags = Reduce(function(x, y) merge(x, y, all=TRUE), list(qc_flags, midqc, sstqc, nbackqc))

# ##################################################
# ########## adding task fmri motion parameters
# ##################################################

# tmp = formattxt('nback_bwroi02.txt', datapath)
# nbackmotion = tmp[,c(fixed_vars, 'tfmri_nback_all_beta_mm', 'tfmri_nback_ab_subthnvols')]
# colnames(nbackmotion)[c(ncol(nbackmotion)-1, ncol(nbackmotion))] = c('tfmri_nback_all_b_meanmotion', 'tfmri_nback_all_b_subthreshnvols')

# tmp = formattxt('mrisst02.txt', datapath)
# sstmotion = tmp[,c(fixed_vars, 'tfmri_sa_beta_mm', 'tfmri_sa_beta_stvols')]
# colnames(sstmotion)[c(ncol(sstmotion)-1, ncol(sstmotion))] = c('tfmri_sst_all_b_meanmotion', 'tfmri_sst_all_b_subthreshnvols')

# tmp = formattxt('midaparc03.txt', datapath) ## N = 14650
# midmotion = tmp[,c(fixed_vars, 'tfmri_mid_all_b_meanmotion', 'tfmri_mid_all_b_subthreshnvols')]

# qc_flags = Reduce(function(x, y) merge(x, y, all=TRUE), list(qc_flags, nbackmotion,sstmotion, midmotion))

# # qc_flags = qc_flags[qc_flags$eventname=='baseline_year_1_arm_1',]

# #####################
# # x = qc_flags
# # tt = x$src_subject_id[duplicated(x$src_subject_id)]
# # tmp = {}
# # for(i in 1:length(tt)){
# # 	if(sum(duplicated(x$eventname[x$src_subject_id==tt[i]])) == 0){
# # 		ind = 1
# # 	}else{
# # 		print(paste0('duplicate ID:', tt[i]))
# # 		tmp = c(tmp,toString(tt[i]))
# # 	}
# # }

# ####### duplicates
# # dups = c('NDAR_INV2ZA2LC3N', 'NDAR_INV3E0WVH3G', 'NDAR_INVJ9GNXGK5', 'NDAR_INVWE1DE80Z', 'NDAR_INVXN6HMGK8')
# # visitiid_toremove = c('S065_INV2ZA2LC3N_20171215', 'S076_INV3E0WVH3G_20171129', 'S053_INVJ9GNXGK5_20180116', 'S076_INVWE1DE80Z_20170913', 'S076_INVXN6HMGK8_20170913')
# # for(i in 1:length(dups)){
# # 	qc_flags[which(qc_flags$src_subject_id ==dups[i]),] ## match mri_info_visitid with imaging data visitid, remove those that are not in imaging data
# # 	qc_flags = qc_flags[-which(qc_flags$mri_info_visitid ==visitiid_toremove[i]), ]
# # }
# # qc_flags[which(qc_flags$src_subject_id ==dups[i]),] ## match mri_info_visitid with imaging data visitid, remove those that are not in imaging data
# # qc_flags = qc_flags[-which(qc_flags$mri_info_visitid =='S076_INVWE1DE80Z_20170913'), ]

# ################# calculate mean motion for taskres_Rs_roi
# # qc_flags$tfmri_alltask_meanmotion = (qc_flags$tfmri_nback_all_b_meanmotion*qc_flags$tfmri_nback_all_b_subthreshnvols + qc_flags$tfmri_sst_all_b_meanmotion*qc_flags$tfmri_sst_all_b_subthreshnvols + qc_flags$tfmri_mid_all_b_meanmotion*qc_flags$tfmri_mid_all_b_subthreshnvols)/(qc_flags$tfmri_nback_all_b_subthreshnvols + qc_flags$tfmri_sst_all_b_subthreshnvols + qc_flags$tfmri_mid_all_b_subthreshnvols)
# # qc_flags$tfmri_alltask_sumnts = qc_flags$tfmri_nback_all_b_subthreshnvols + qc_flags$tfmri_sst_all_b_subthreshnvols + qc_flags$tfmri_mid_all_b_subthreshnvols

# write.csv(qc_flags, paste0('/home/wez025/matlab/ABCD_release_3.0/Released/Data/', 'qcflags_dof_nts_meanmotion_release4.0.csv'), row.names = FALSE)

# tt = read.csv(paste0('/home/wez025/matlab/ABCD_release_3.0/Released/Data/', 'qcflags_dof_nts_meanmotion_release4.0.csv'))
# )
# tt = tt[tt$eventname=='baseline_year_1_arm_1', ]
# tt[11817,]
# ######################################################################################################
# ######################################################################################################
# ######################################################################################################
# ################################## COMBINING: behavioral data extraction
# ######################################################################################################
# ######################################################################################################
# ######################################################################################################
# fixed_vars = c("src_subject_id",'eventname','interview_date','interview_age','sex')

# tmp = formattxt('abcd_saag01.txt', datapath)
# sag = tmp[,c(fixed_vars, 'sag_miss_school_excuse_p', 'sag_miss_school_unexcuse_p', 'sag_excused_absence_p', 'sag_unexcuse_absence_p', 'sag_grade_type', 'sag_iep_current_p', 'sag_iep_ever_p')]
# #### sag_grade_type: -1=not applicable; 777=refuse to answer; 1=A+, 12=F


# ##### CBCL
# vars = c('cbcl_scr_syn_aggressive_r',
# 'cbcl_scr_syn_anxdep_r',
# 'cbcl_scr_syn_attention_r',
# 'cbcl_scr_syn_rulebreak_r',
# 'cbcl_scr_syn_social_r',
# 'cbcl_scr_syn_somatic_r',
# 'cbcl_scr_syn_thought_r',
# 'cbcl_scr_syn_withdep_r',
# 'cbcl_scr_syn_external_r',
# 'cbcl_scr_syn_internal_r',
# 'cbcl_scr_syn_totprob_r'
# )
# cbcl = formattxt('abcd_cbcls01.txt', datapath)
# cbcl = cbcl[, c(fixed_vars, vars)]

# ###### bpm & poa
# bpmpoa = formattxt('abcd_yssbpm01.txt', datapath)
# vars = c('bpm_y_ss_attention_mean',
# 'bpm_y_ss_external_mean',
# 'bpm_y_ss_internal_mean',
# 'bpm_y_ss_totalprob_mean',
# 'poa_y_ss_sum'
# )
# bpmpoa = bpmpoa[,c(fixed_vars, vars)]

# ##### upps & bisbas
# vars = c('upps_y_ss_lack_of_perseverance',
# 'upps_y_ss_lack_of_planning',
# 'upps_y_ss_negative_urgency',
# 'upps_y_ss_positive_urgency',
# 'upps_y_ss_sensation_seeking',
# 'bis_y_ss_bas_drive',
# 'bis_y_ss_bas_fs',
# 'bis_y_ss_bas_rr',
# 'bis_y_ss_bis_sum'
# )
# uppsbisbas = formattxt('abcd_mhy02.txt', datapath)
# uppsbisbas = uppsbisbas[,c(fixed_vars, vars)]

# ##### NIH TB
# vars = c(
# 'nihtbx_flanker_uncorrected',
# 'nihtbx_list_uncorrected',
# 'nihtbx_cardsort_uncorrected',
# 'nihtbx_reading_uncorrected',
# 'nihtbx_pattern_uncorrected',
# 'nihtbx_picture_uncorrected',
# 'nihtbx_picvocab_uncorrected',
# 'nihtbx_cryst_uncorrected',
# 'nihtbx_fluidcomp_uncorrected',
# 'nihtbx_totalcomp_uncorrected')

# nihtb = formattxt('abcd_tbss01.txt', datapath)
# nihtb = nihtb[,c(fixed_vars, vars)]

# ###### little man task
# ###### efficiency = lmt_scr_perc_correct/lmt_scr_rt_correct
# vars = c('lmt_scr_perc_correct',
# 'lmt_scr_rt_correct')
# lmt = formattxt('lmtp201.txt',datapath)
# lmt = lmt[,c(fixed_vars, vars)] 

# ###### matrix reasoning & RAVLT
# vars = c('pea_wiscv_trs',
# 'pea_ravlt_sd_trial_i_tc' ,'pea_ravlt_sd_trial_ii_tc','pea_ravlt_sd_trial_iii_tc','pea_ravlt_sd_trial_iv_tc','pea_ravlt_sd_trial_v_tc',
# 'pea_ravlt_ld_trial_vii_tc')
# ps = formattxt('abcd_ps01.txt', datapath)
# ps = ps[,c(fixed_vars, vars)]

# ###### emotional stroop
# vars = c('strp_scr_mnrt_congr',
# 'strp_scr_mnrt_incongr', 
# 'strp_scr_mnrt_all',
# 'strp_scr_acc_congr',
# 'strp_scr_acc_incongr',
# 'strp_scr_acc_all'
# )
# strp = formattxt('abcd_yest01.txt', datapath)
# strp = strp[, c(fixed_vars, vars)]

# ##### delayed discounting task: only data from half of the sample are available

# vars = c('ddis_scr_val_indif_point_6h',
# 'ddis_scr_val_indif_pnt_1da',
# 'ddis_scr_val_indif_pnt_1week',
# 'ddis_scr_val_indif_pnt_1mth',
# 'ddis_scr_val_indif_pnt_3mth',
# 'ddis_scr_val_indif_pnt_1yr',
# 'ddis_scr_val_indif_pnt_5yr'
# )
# dds = formattxt('abcd_yddss01.txt', datapath)
# dds = dds[,c(fixed_vars, vars)]

# vars = c('tfmri_sst_all_beh_total_mssrt')
# ssrt = formattxt('abcd_sst02.txt', datapath)
# ssrt = ssrt[, c(fixed_vars, vars)]
# names(ssrt)[dim(ssrt)[2]] = 'ssrt'

# behVars = Reduce(function(x, y) merge(x, y, all=TRUE), list(dds,strp, ps, lmt, nihtb, cbcl, uppsbisbas, bpmpoa, ssrt))

# x = behVars
# tt = x$src_subject_id[duplicated(x$src_subject_id)]
# tmp = {}
# for(i in 1:length(tt)){
# 	if(sum(duplicated(x$eventname[x$src_subject_id==tt[i]])) == 0){
# 		ind = 1
# 	}else{
# 		print(paste0('duplicate ID:', tt[i]))
# 		tmp = c(tmp,toString(tt[i]))
# 	}
# }
# write.csv(behVars, paste0('/home/wez025/matlab/ABCD_release_3.0/Released/Data/', 'behVars_release4.0.csv'), row.names = FALSE)