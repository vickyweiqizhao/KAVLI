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

# dat = read.csv('/home/wez025/matlab/ABCD_release_3.0/Released/Data/qcflags_dof_nts_meanmotion_release4.0.csv')
# dat = dat[,c("src_subject_id", 'interview_date', 'interview_age', 'eventname', 'mri_info_visitid', 'imgincl_nback_include')]
# # dat = dat[which(dat$imgincl_nback_include == 1),]
# dat = dat[which(dat$eventname == 'baseline_year_1_arm_1'),]

# write.csv(dat, paste0('/home/wez025/matlab/ABCD_code/KAVLI/data/', 'subjects_passQC_4.0.csv'), row.names = FALSE)

######################################################################################################
####################### nback ROI data
######################################################################################################
fixed_vars = c('src_subject_id','eventname','interview_date', 'interview_age', 'sex')

tmp = formattxt('nbackr201.txt', datapath) ## second run nback
tmp = tmp[which(tmp$eventname=='baseline_year_1_arm_1'),]

nback_ROI = tmp[,c(fixed_vars, names(tmp)[which(names(tmp)=='tfmri_nback_r2_121'):which(names(tmp)=='tfmri_nback_r2_150')], names(tmp)[which(names(tmp)=='tfmri_nback_r2_543'):which(names(tmp)=='tfmri_nback_r2_610')])]
nback_ROI$eventname = 'baseline'

nback_ROI = nback_ROI[complete.cases(nback_ROI),]
nback_ROI$src_subject_id = substring(nback_ROI$src_subject_id, 6)

outputdir = '/home/wez025/matlab/KAVLI_data/timeseries_by_ROI/'
write.csv(nback_ROI, paste0(outputdir, 'nback_ROI.csv'), row.names = FALSE)
