library(R.matlab)
library(readxl)
library(dplyr)
library(tidyr)

####################### paths and importing functions
datapath = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/tabulated/released/'
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
# dat = dat[which(dat$imgincl_nback_include == 1),]
dat = dat[which(dat$eventname == 'baseline_year_1_arm_1'),]

write.csv(dat, paste0('/home/wez025/matlab/ABCD_code/KAVLI/data/', 'subjects_passQC_4.0.csv'), row.names = FALSE)

######################################################################################################
####################### sMRI data
######################################################################################################
fixed_vars = c('src_subject_id','eventname','interview_date', 'interview_age', 'sex')

tmp = formattxt('abcd_smrip10201.txt', datapath) ## N = 17485
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
thickness = names(tmp)[grep("smri_thick_cdk", names(tmp))]
area = names(tmp)[grep("smri_area_cdk", names(tmp))]
sulc = names(tmp)[grep("smri_sulc_cdk", names(tmp))]
vol = names(tmp)[grep("smri_vol_cdk", names(tmp))]
smri = tmp[,c(fixed_vars, thickness, area, sulc, vol)]
means = names(smri)[grep("mean", names(smri))]
total = names(smri)[grep("total", names(smri))]
smri = smri[,-c(grep("mean", names(smri)), grep("total", names(smri)))]

tmp = formattxt('abcd_smrip20201.txt', datapath) ## T1 intensity measures
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
t1w = tmp[,non_means]
t1w = t1w[, setdiff(names(t1w), c("collection_title", "study_cohort_name"))]

tmp = formattxt('abcd_smrip30201.txt', datapath) ## T2 intensity measures
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
t2w = tmp[,non_means]
t2w = t2w[, setdiff(names(t2w), c("collection_title", "study_cohort_name", 'abcd_smrip30201_id', 'subjectkey'))]

structural = Reduce(function(x, y) merge(x, y, all=TRUE), list(smri, t1w,t2w))

####### DTI data

tmp = formattxt('abcd_dti_p101.txt', datapath) ## T2 intensity measures
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
dti1 = tmp[,non_means]
dti1 = dti1[, setdiff(names(dti1), c("collection_title", "study_cohort_name", 'abcd_dti_p101_id', 'subjectkey'))] ### contains visitID!!
names(dti1)[which(names(dti1)=='dmri_dti_visitid')] = 'mri_visitID'

tmp = formattxt('abcd_dti_p201.txt', datapath) ## T2 intensity measures
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
dti2 = tmp[,non_means]
dti2 = dti2[, setdiff(names(dti2), c("collection_title", "study_cohort_name", 'abcd_dti_p201_id', 'subjectkey'))] ### contains visitID!!

dti = Reduce(function(x, y) merge(x, y, all=TRUE), list(dti1, dti2))

####### RSI data

tmp = formattxt('abcd_drsip101.txt', datapath) ## RSI
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
rsi1 = tmp[,non_means]
rsi1 = rsi1[, setdiff(names(rsi1), c("collection_title", "study_cohort_name", 'abcd_drsip101_id', 'subjectkey'))] ### contains visitID!!
names(rsi1)[which(names(rsi1)=='dmri_rsi_visitid')] = 'mri_visitID'

tmp = formattxt('abcd_drsip201.txt', datapath) ## RSI
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
rnd = tmp[,non_means]
rnd = rnd[, setdiff(names(rnd), c("collection_title", "study_cohort_name", 'abcd_drsip201_id', 'subjectkey'))] ### contains visitID!!

tmp = formattxt('abcd_drsip301.txt', datapath) ## RSI
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
rnt = tmp[,non_means]
rnt = rnt[, setdiff(names(rnt), c("collection_title", "study_cohort_name", 'abcd_drsip301_id', 'subjectkey'))] ### contains visitID!!

tmp = formattxt('abcd_drsip401.txt', datapath) ## RSI
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
hni = tmp[,non_means]
hni = hni[, setdiff(names(hni), c("collection_title", "study_cohort_name", 'abcd_drsip401_id', 'subjectkey'))] ### contains visitID!!

tmp = formattxt('abcd_drsip501.txt', datapath) ## RSI
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
hnd = tmp[,non_means]
hnd = hnd[, setdiff(names(hnd), c("collection_title", "study_cohort_name", 'abcd_drsip501_id', 'subjectkey'))] ### contains visitID!!

tmp = formattxt('abcd_drsip601.txt', datapath) ## RSI
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
hnt = tmp[,non_means]
hnt = hnt[, setdiff(names(hnt), c("collection_title", "study_cohort_name", 'abcd_drsip601_id', 'subjectkey'))] ### contains visitID!!

tmp = formattxt('abcd_drsip701.txt', datapath) ## RSI
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]
non_means = setdiff(names(tmp), names(tmp)[grep("mean", names(tmp))])
fni = tmp[,non_means]
fni = fni[, setdiff(names(fni), c("collection_title", "study_cohort_name", 'abcd_drsip701_id', 'subjectkey'))] ### contains visitID!!

rsi = Reduce(function(x, y) merge(x, y, all=TRUE), list(fni, hni, rnt, rsi1, rnd))

###### concatenate everything

qcflag = read.csv('/home/wez025/matlab/ABCD_release_3.0/Released/Data/qcflags_dof_nts_meanmotion_release4.0.csv')
qcflag = qcflag[which(qcflag$eventname == 'baseline_year_1_arm_1'),]
qcflag$rsi_dti_smri_pass = qcflag$imgincl_t1w_include + qcflag$imgincl_t2w_include + qcflag$imgincl_dmri_include + qcflag$fsqc_qc
qcpass = qcflag[,c(fixed_vars, "rsi_dti_smri_pass")]

# structuralpass = Reduce(function(x, y) merge(x, y, all=TRUE), list(qcflag, structural))
# structuralpass = structuralpass[intersect(which(structuralpass$imgincl_t1w_include == 1), which(structuralpass$imgincl_t2w_include == 1)), ]
# structuralpass

data = Reduce(function(x, y) merge(x, y, all=TRUE), list(qcpass, rsi, dti, structural))
data = data[,-which(names(data)=="subjectkey")]
data = data[data$rsi_dti_smri_pass==4, ]
nancols = names(data[1,is.na(data[1,])])
data = data[,setdiff(names(data), nancols)]
data = data[complete.cases(data),]
data = data[,-which(names(data)=="rsi_dti_smri_pass")]
names(data)[which(names(data)=='mri_visitID')] = 'mri_info_visitid'

structure_cols = c('mri_info_visitid', names(data)[grep('smri_', names(data))])
structural_dat = data[,structure_cols]

rsi_cols = c('mri_info_visitid', names(data)[grep('dmri_rsi', names(data))])
rsi_dat = data[,rsi_cols]

dti_cols = c('mri_info_visitid', names(data)[grep('dmri_dti', names(data))])
dti_dat = data[,dti_cols]

demographics = data[, 1:6]

outputdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/'
write.csv(structural_dat, paste0(outputdir, 'structure_data.csv'))
write.csv(dti_dat, paste0(outputdir, 'dti_data.csv'))
write.csv(rsi_dat, paste0(outputdir, 'rsi_data.csv'))
write.csv(demographics, paste0(outputdir, 'demographics_data.csv'))