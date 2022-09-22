#########Preprocessing for GTEx EPIC Data#########
######Run in R 4.0.3######
#setwd("U:/")
.libPaths(new="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/R_4.0_methylation")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
##BiocManager::install()
#library(BiocManager)
#BiocManager::install("minfi")
#BiocManager::install("ChAMP")
#BiocManager::install("sva")
#BiocManager::install("limma")
#BiocManager::install("missMethyl")
#BiocManager::install("wateRmelon")
#BiocManager::install("WGCNA")
#BiocManager::install("DMRcate")
#BiocManager::install("variancePartition")
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b3.hg19"

######Load libraries needed for preprocessing######
library("ChAMP")
library("minfi")
library("IlluminaHumanMethylationEPICanno.ilm10b3.hg19")
library("stringr")
library("wateRmelon")

######EPIC Manifest file from Illumina######
manifest <- read.csv("G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/MethylationEPIC_v-1-0_B4.csv", header=T, stringsAsFactors=F)
dim(manifest)

######Platemap for EPIC array experiment######
###Plate number and position for each sample.
meth_platemap <- read.csv("G:/Projects/GTEx/Telomere/eGTExDNA_Pierce_Jan'18.csv", stringsAsFactors=FALSE)
dim(meth_platemap)
head(meth_platemap)

###Rename and reformat GTEx IDs for platemap.
meth_platemap$CollaboratorSampleID <- meth_platemap$Collaborator.Sample.ID
table(nchar(meth_platemap$CollaboratorSampleID)) ###all 14 or 15 characters, so add SM 
meth_platemap$CollaboratorSampleID[nchar(meth_platemap$CollaboratorSampleID)==15] <- paste(meth_platemap$CollaboratorSampleID[nchar(meth_platemap$CollaboratorSampleID)==15], "-SM", sep="")
meth_platemap$CollaboratorSampleID[nchar(meth_platemap$CollaboratorSampleID)==14] <- paste(meth_platemap$CollaboratorSampleID[nchar(meth_platemap$CollaboratorSampleID)==14], "-SM", sep="")

length(unique(meth_platemap$CollaboratorSampleID)) ###1000 unique samples no duplicates 
length(unique(meth_platemap$Tissue.Site.Detail))  ###9 tissues 
length(unique(meth_platemap$Collaborator.Participant.ID)) ##428 donors
table(meth_platemap$Tissue.Site.Detail)

######GTEX covariate file######
gtex_cov <- read.csv("G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/gtex_covariates_all_samples.csv", header=TRUE, stringsAsFactors=FALSE)
dim(gtex_cov)
head(gtex_cov)

######Merge methylation platemap and GTEX covariate file######
meth_platemap_cov <- merge(meth_platemap, gtex_cov, by="CollaboratorSampleID", all.x=TRUE)
dim(meth_platemap_cov) ###1000 
summary(meth_platemap_cov) ###No missing donor information.  

######Preprocessing using ChAMP pipeline######
###Identify base directory that holds .idat files. 
#baseDir <-file.path("G:/Projects/GTEx/GTEX_EPIC_Meth_Data/")
#list.files(baseDir)

###Load .idat files. Only filtered poorly performing probes and samples at this step.  
#raw <- champ.load(directory=baseDir, method="minfi", arraytype="EPIC", filterDetP=TRUE,
#ProbeCutoff=0, SampleCutoff=0.05, detPcut=0.01, filterBeads=TRUE, beadCutoff=0.05,filterNoCG=TRUE,
#filterSNPs=FALSE, filterMultiHit=FALSE, filterXY=FALSE)

###Save ChAMP data object. Contains raw methylation data, pd sample sheet, and other information.
#save(raw, file="U:/GTEX_methylation/CHAMP_raw_8-8-2019.RData")

######Load ChAMP data object (named raw)###### 
load(file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/CHAMP_raw_8-8-2019.RData")
raw$pd <- raw$pd[order(raw$pd$Sample_Name), ]
raw$beta <- raw$beta[ ,order(colnames(raw$beta))]

###Update pd sample sheet (contains covariates related to methylation experiment).
raw$pd$good <- 1 ###samples that passed filtering.
###Extracted information related to Container and position.  
subIDs <- as.data.frame(str_split_fixed(as.character(raw$pd$Sample_Name), "_", 3))
colnames(subIDs) <- c("id1", "id2", "id3")
raw$pd$Container <- as.character(subIDs$id2)
raw$pd$Position <- as.character(subIDs$id3)

######Merge platemap+gtex covariate dataset with pd sample sheet####### 
meth_data <- merge(raw$pd, meth_platemap_cov, by=c("Container", "Position"))
dim(meth_data) ###997 x 44 
meth_platemap_cov$CollaboratorSampleID[!(meth_platemap_cov$CollaboratorSampleID %in% meth_data$CollaboratorSampleID)]
###3 samples were poorly performing:GTEX-12WSI-0826-SM, GTEX-1PBJI-0826-SM, GTEX-UTHO-0003-SM

#####Filter samples that are potential sex mismatches###### 
raw_gs <- mapToGenome(raw$rgSet)  
raw_cn <- getCN(raw_gs) ###extracts copy number estimates from probes 

x <- manifest$Name[which(manifest$CHR=="X")]
length(x)
y <- manifest$Name[which(manifest$CHR=="Y")]
length(y)
x_keep <- rownames(raw$beta)[rownames(raw$beta) %in% x]
length(x_keep) ##17197...enough to use for x. 
y_keep <- rownames(raw$beta)[rownames(raw$beta) %in% y]
length(y_keep) ##84 too few to discern males from females. Using full set of y probes.  

raw_cn_x <- raw_cn[rownames(raw_cn) %in% x_keep, ]
raw_cn_x_medians <- colMedians(raw_cn_x, na.rm=TRUE)
raw_cn_y <- raw_cn[rownames(raw_cn) %in% y, ]
raw_cn_y_medians <- colMedians(raw_cn_y, na.rm=TRUE)
differences <- raw_cn_y_medians - raw_cn_x_medians
hist(differences)

sex_prediction <- DataFrame(xMed=raw_cn_x_medians, yMed=raw_cn_y_medians, difference=differences)
rownames(sex_prediction) <- colnames(raw_cn) 
sex_prediction$Sample_Name <- rownames(sex_prediction)

meth_data <- merge(meth_data, sex_prediction, by="Sample_Name")
dim(meth_data)
meth_data$predictedSex[meth_data$difference <= -2] <- "F"
meth_data$predictedSex[meth_data$difference > -2]  <- "M"

meth_data$sex_mismatch <- 0
meth_data$sex_mismatch[meth_data$predictedSex=="F" & meth_data$SEX==1] <- 1
meth_data$sex_mismatch[meth_data$predictedSex=="M" & meth_data$SEX==2] <- 1
table(meth_data$sex_mismatch)
mismatch_sex_samples <- meth_data$CollaboratorSampleID[meth_data$sex_mismatch==1]
###Will remove 6 samples due to potential sex mismatch. Just to be on conservative side:GTEX-14PKV-1226-SM,  GTEX-14PN4-0626-SM, GTEX-14PII-0626-SM, GTEX-15DCE-0426-SM,  GTEX-13N11-1726-SM, GTEX-1212Z-1226-SM 
meth_data <- meth_data[!(meth_data$CollaboratorSampleID %in% mismatch_sex_samples), ]
dim(meth_data)  ###991 x 49 

######Filter samples based on mismatched genotype###### 
beta_minfi_snp <- getSnpBeta(raw$rgSet)
dim(beta_minfi_snp)
beta_minfi_snp[1:3, 1:3]
snps  <- matrix(unlist(genme(beta_minfi_snp, peaks=c(0.2, 0.5, 0.8))), nrow=59, ncol=997, byrow=FALSE)
dim(snps)
snps[1:3, 1:3]
beta_minfi_snp[1:3, 1:3]

geno_snps <- data.frame(snps)
colnames(geno_snps) <- colnames(beta_minfi_snp)
rownames(geno_snps) <- rownames(beta_minfi_snp)
geno_snps <- data.frame(t(geno_snps))
geno_snps$meth_sample <- rownames(geno_snps)
geno_snps[geno_snps==1] <- 0 
geno_snps[geno_snps==2] <- 1
geno_snps[geno_snps==3] <- 2
head(geno_snps)

meth_data$meth_sample <- meth_data$Sample_Name
genotype_review_ids <- meth_data[ ,c("meth_sample", "CollaboratorSampleID", "CollaboratorParticipantID", "Tissue.Site.Detail")]
geno_snps_ids <- merge(genotype_review_ids, geno_snps, by="meth_sample")
head(geno_snps_ids)
dim(geno_snps_ids)
geno_snps_ids <- geno_snps_ids[order(geno_snps_ids$CollaboratorSampleID), ]
head(geno_snps_ids)
#write.csv(geno_snps_ids, "U:/GTEX_methylation/GTEX_samples_genotypes_59SNPs_8-19-2019.csv")

###From Lin
#Obs	sample	ID	Tissue	sum_diff
#170	GTEX-13N1W-2126-SM	GTEX-13N1W	Kidney-Cortex	5
#171	GTEX-13N1W-2326-SM	GTEX-13N1W	Colon-Transvers	5
#172	GTEX-13N1W-2526-SM	GTEX-13N1W	Prostate	5
#541	GTEX-1MCC2-0726-SM	GTEX-1MCC2	Lung	10
#830	GTEX-YJ89-0526-SM	GTEX-YJ89	Muscle-Skeletal	33
####We decied to remove only GTEX-YJ89-0526-SM 
bad_genotype_sample <- "GTEX-YJ89-0526-SM"
meth_data <- meth_data[!(meth_data$CollaboratorSampleID %in% bad_genotype_sample), ]
dim(meth_data)  ###990 x 50 
 
######Filtering probes based on Pidsley et al 2016 documentation######
meth_raw <- raw$beta
dim(meth_raw)
###Removing cross-reactive probes###
cross <- read.csv("G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/crossreactives_pidsley.csv", header=T) 
dim(cross) ##43254 probes
meth_raw_cross <- meth_raw[!(rownames(meth_raw) %in% cross$X), ]
dim(meth_raw_cross) ###778033 probes. 

###SNPs within single base extension SBE or in CpG###
cpg_snp <- read.csv("G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/snp_cpg_pidsley.csv", header=T)
dim(cpg_snp) ##12510 probes
meth_raw_cross_snps <- meth_raw_cross[!(rownames(meth_raw_cross) %in% cpg_snp$PROBE), ]
dim(meth_raw_cross_snps) ###770541.  

sbe_snp <- read.csv("G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/snp_sbe_pidsley.csv", header=T)
dim(sbe_snp) ##414 probes
meth_raw_cross_snps <- meth_raw_cross_snps[!(rownames(meth_raw_cross_snps) %in% sbe_snp$PROBE), ]
dim(meth_raw_cross_snps) ###770325

###Remove RS probes. ChAMP already filters these probes from beta dataframe.### 
rs <- manifest$Name[grep("rs", manifest$Name)]
length(rs) ##59
meth_raw_cross_snps_rs <- meth_raw_cross_snps[!(rownames(meth_raw_cross_snps) %in% rs), ]
dim(meth_raw_cross_snps_rs)##770325.

###Remove X and Y probes### 
meth_raw_cross_snps_rs_X <- meth_raw_cross_snps_rs[!(rownames(meth_raw_cross_snps_rs) %in% x), ]
dim(meth_raw_cross_snps_rs_X) ##754313
meth_raw_cross_snps_rs_XY <- meth_raw_cross_snps_rs_X[!(rownames(meth_raw_cross_snps_rs_X) %in% y), ]
dim(meth_raw_cross_snps_rs_XY) ###754288

###Remove remaining probes that are poor performing according to Illumina (update from B4 to B3)###  
cpgs_manifest_good <- rownames(meth_raw_cross_snps_rs_XY)[rownames(meth_raw_cross_snps_rs_XY) %in% manifest$Name]
length(cpgs_manifest_good) ###754119
meth_raw_cross_snps_rs_XY_good <- meth_raw_cross_snps_rs_XY[rownames(meth_raw_cross_snps_rs_XY) %in% cpgs_manifest_good, ]
dim(meth_raw_cross_snps_rs_XY_good) ###754119


######Preprocessing methylation dataset using background correction and BMIQ normalization######
###covariate file
cov <- data.frame(meth_data)
dim(cov)
cov <- cov[order(cov$Sample_Name), ]
head(cov)
final_ids <- cov$Sample_Name
length(final_ids)
final_ids[1:6]

###BMIQ only### 
meth_final <- meth_raw_cross_snps_rs_XY_good[ ,colnames(meth_raw_cross_snps_rs_XY_good) %in% final_ids]
meth_final <- meth_final[ ,order(colnames(meth_final))]
dim(meth_final) 
set.seed(198833) #not needed anymore
meth_final_BMIQ <- champ.norm(beta=meth_final, method="BMIQ", plotBMIQ=FALSE, arraytype="EPIC",cores=1)
colnames(meth_final_BMIQ) <- cov$CollaboratorSampleID
dim(meth_final_BMIQ)
save(meth_final_BMIQ,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/raw_beta_final_BMIQ_990_2-4-2021.RData")

###raw datafile to output 
colnames(meth_final) <- cov$CollaboratorSampleID
dim(meth_final)
save(meth_final,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/raw_beta_final_990_2-4-2021.RData")

###Preprocess using Noob and BMIQ### 
noob <- preprocessNoob(raw$rgSet, dyeCorr = TRUE, verbose = FALSE, dyeMethod="single")
noob_beta <- getBeta(noob)
noob_beta_final <- noob_beta[rownames(noob_beta) %in% rownames(meth_final), ]
dim(noob_beta_final)
noob_beta_final <- noob_beta_final[,colnames(noob_beta_final) %in% final_ids]
noob_beta_final <- noob_beta_final[ ,order(colnames(noob_beta_final))]
set.seed(198833) #not needed anymore
noob_beta_final_BMIQ <- champ.norm(beta=noob_beta_final, method="BMIQ", plotBMIQ=FALSE, arraytype="EPIC",cores=1)
colnames(noob_beta_final_BMIQ) <- cov$CollaboratorSampleID
save(noob_beta_final_BMIQ,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/noob_beta_final_BMIQ_990_2-4-2021.RData")

###Preprocess using Illumina background correction and BMIQ### 
bc <- preprocessIllumina(raw$rgSet, bg.correct = TRUE, normalize="no")
bc_beta <- getBeta(bc)
bc_beta_final <- bc_beta[rownames(bc_beta) %in% rownames(meth_final), ]
dim(bc_beta_final)
bc_beta_final <- bc_beta_final[,colnames(bc_beta_final) %in% final_ids]
bc_beta_final <- bc_beta_final[,order(colnames(bc_beta_final))]
set.seed(198833) #not needed anymore
bc_beta_final_BMIQ <- champ.norm(beta=bc_beta_final, method="BMIQ", plotBMIQ=FALSE, arraytype="EPIC",cores=1)
colnames(bc_beta_final_BMIQ) <- cov$CollaboratorSampleID
save(bc_beta_final_BMIQ,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/GTEX_methylation_preprocessing/bc_beta_final_BMIQ_990_2-4-2021.RData")


######Tissue specific datasets###### 
keep_covariates <- c("CollaboratorSampleID", "CollaboratorParticipantID", "TSD", "Sample_Name", "Container", "Position", "Sample_Plate", "Array","Slide", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY",  "sample_RINmean", "TQImean", "raceRE", "ALS", "Autolysis.Score")
cov_final_selected <- cov[ ,colnames(cov) %in% keep_covariates]
dim(cov_final_selected)

###Function to generate datasets for each tissue### 
tissue_specific_methylation <- function(tissue, methylation_data)
{
ids_tissue <- cov$CollaboratorSampleID[cov$TSD==tissue]
methylation_tissue <- methylation_data[ ,colnames(methylation_data) %in% ids_tissue]
methylation_tissue <- methylation_tissue[,order(colnames(methylation_tissue))]
methylation_tissue <- methylation_tissue[order(rownames(methylation_tissue)), ]
return(methylation_tissue)
}

###Whole Blood### 
raw_final_BMIQ_wb <- tissue_specific_methylation("Whole Blood", meth_final_BMIQ)
dim(raw_final_BMIQ_wb)
bc_final_BMIQ_wb <- tissue_specific_methylation("Whole Blood", bc_beta_final_BMIQ)
dim(bc_final_BMIQ_wb)
noob_final_BMIQ_wb <- tissue_specific_methylation("Whole Blood", noob_beta_final_BMIQ)
dim(noob_final_BMIQ_wb)
save(raw_final_BMIQ_wb,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/raw_beta_bmiq_2-6-2021/raw_final_BMIQ_wb_2-6-2021.RData")
save(bc_final_BMIQ_wb,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/bc_beta_bmiq_2-6-2021/bc_final_BMIQ_wb_2-6-2021.RData")
save(noob_final_BMIQ_wb,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/noob_beta_bmiq_2-6-2021/noob_final_BMIQ_wb_2-6-2021.RData")

cov_wb <- cov_final_selected[which(cov_final_selected$TSD=="Whole Blood"), ]
cov_wb <- cov_wb[order(cov_wb$CollaboratorSampleID), ]
write.csv(cov_wb, "G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/covariates_tissue_2-6-2021/wb_covariates_2-6-2021.csv", row.names=FALSE)

###Lung###
raw_final_BMIQ_lung <- tissue_specific_methylation("Lung", meth_final_BMIQ)
dim(raw_final_BMIQ_lung)
bc_final_BMIQ_lung <- tissue_specific_methylation("Lung", bc_beta_final_BMIQ)
dim(bc_final_BMIQ_lung)
noob_final_BMIQ_lung <- tissue_specific_methylation("Lung", noob_beta_final_BMIQ)
dim(noob_final_BMIQ_lung)
save(raw_final_BMIQ_lung,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/raw_beta_bmiq_2-6-2021/raw_final_BMIQ_lung_2-6-2021.RData")
save(bc_final_BMIQ_lung,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/bc_beta_bmiq_2-6-2021/bc_final_BMIQ_lung_2-6-2021.RData")
save(noob_final_BMIQ_lung,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/noob_beta_bmiq_2-6-2021/noob_final_BMIQ_lung_2-6-2021.RData")

cov_lung <- cov_final_selected[which(cov_final_selected$TSD=="Lung"), ]
cov_lung <- cov_lung[order(cov_lung$CollaboratorSampleID), ]
write.csv(cov_lung, "G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/covariates_tissue_2-6-2021/lung_covariates_2-6-2021.csv", row.names=FALSE)

###Prostate###
raw_final_BMIQ_prostate <- tissue_specific_methylation("Prostate", meth_final_BMIQ)
dim(raw_final_BMIQ_prostate)
bc_final_BMIQ_prostate <- tissue_specific_methylation("Prostate", bc_beta_final_BMIQ)
dim(bc_final_BMIQ_prostate)
noob_final_BMIQ_prostate <- tissue_specific_methylation("Prostate", noob_beta_final_BMIQ)
dim(noob_final_BMIQ_prostate)
save(raw_final_BMIQ_prostate,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/raw_beta_bmiq_2-6-2021/raw_final_BMIQ_prostate_2-6-2021.RData")
save(bc_final_BMIQ_prostate,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/bc_beta_bmiq_2-6-2021/bc_final_BMIQ_prostate_2-6-2021.RData")
save(noob_final_BMIQ_prostate,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/noob_beta_bmiq_2-6-2021/noob_final_BMIQ_prostate_2-6-2021.RData")

cov_prostate <- cov_final_selected[which(cov_final_selected$TSD=="Prostate"), ]
cov_prostate <- cov_prostate[order(cov_prostate$CollaboratorSampleID), ]
write.csv(cov_prostate, "G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/covariates_tissue_2-6-2021/prostate_covariates_2-6-2021.csv", row.names=FALSE)

###Testis###
raw_final_BMIQ_testis <- tissue_specific_methylation("Testis", meth_final_BMIQ)
dim(raw_final_BMIQ_testis)
bc_final_BMIQ_testis <- tissue_specific_methylation("Testis", bc_beta_final_BMIQ)
dim(bc_final_BMIQ_testis)
noob_final_BMIQ_testis <- tissue_specific_methylation("Testis", noob_beta_final_BMIQ)
dim(noob_final_BMIQ_testis)
save(raw_final_BMIQ_testis,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/raw_beta_bmiq_2-6-2021/raw_final_BMIQ_testis_2-6-2021.RData")
save(bc_final_BMIQ_testis,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/bc_beta_bmiq_2-6-2021/bc_final_BMIQ_testis_2-6-2021.RData")
save(noob_final_BMIQ_testis,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/noob_beta_bmiq_2-6-2021/noob_final_BMIQ_testis_2-6-2021.RData")

cov_testis <- cov_final_selected[which(cov_final_selected$TSD=="Testis"), ]
cov_testis <- cov_testis[order(cov_testis$CollaboratorSampleID), ]
write.csv(cov_testis, "G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/covariates_tissue_2-6-2021/testis_covariates_2-6-2021.csv", row.names=FALSE)

###Ovary### 
raw_final_BMIQ_ovary <- tissue_specific_methylation("Ovary", meth_final_BMIQ)
dim(raw_final_BMIQ_ovary)
bc_final_BMIQ_ovary <- tissue_specific_methylation("Ovary", bc_beta_final_BMIQ)
dim(bc_final_BMIQ_ovary)
noob_final_BMIQ_ovary <- tissue_specific_methylation("Ovary", noob_beta_final_BMIQ)
dim(noob_final_BMIQ_ovary)
save(raw_final_BMIQ_ovary,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/raw_beta_bmiq_2-6-2021/raw_final_BMIQ_ovary_2-6-2021.RData")
save(bc_final_BMIQ_ovary,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/bc_beta_bmiq_2-6-2021/bc_final_BMIQ_ovary_2-6-2021.RData")
save(noob_final_BMIQ_ovary,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/noob_beta_bmiq_2-6-2021/noob_final_BMIQ_ovary_2-6-2021.RData")

cov_ovary <- cov_final_selected[which(cov_final_selected$TSD=="Ovary"), ]
cov_ovary <- cov_ovary[order(cov_ovary$CollaboratorSampleID), ]
write.csv(cov_ovary, "G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/covariates_tissue_2-6-2021/ovary_covariates_2-6-2021.csv", row.names=FALSE)

###Colon###
raw_final_BMIQ_colon <- tissue_specific_methylation("Colon - Transverse", meth_final_BMIQ)
dim(raw_final_BMIQ_colon)
bc_final_BMIQ_colon <- tissue_specific_methylation("Colon - Transverse", bc_beta_final_BMIQ)
dim(bc_final_BMIQ_colon)
noob_final_BMIQ_colon <- tissue_specific_methylation("Colon - Transverse", noob_beta_final_BMIQ)
dim(noob_final_BMIQ_colon)
save(raw_final_BMIQ_colon,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/raw_beta_bmiq_2-6-2021/raw_final_BMIQ_colon_2-6-2021.RData")
save(bc_final_BMIQ_colon,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/bc_beta_bmiq_2-6-2021/bc_final_BMIQ_colon_2-6-2021.RData")
save(noob_final_BMIQ_colon,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/noob_beta_bmiq_2-6-2021/noob_final_BMIQ_colon_2-6-2021.RData")

cov_colon <- cov_final_selected[which(cov_final_selected$TSD=="Colon - Transverse"), ]
cov_colon <- cov_colon[order(cov_colon$CollaboratorSampleID), ]
write.csv(cov_colon, "G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/covariates_tissue_2-6-2021/colon_covariates_2-6-2021.csv", row.names=FALSE)

###Breast###
raw_final_BMIQ_breast <- tissue_specific_methylation("Breast - Mammary Tissue", meth_final_BMIQ)
dim(raw_final_BMIQ_breast)
bc_final_BMIQ_breast <- tissue_specific_methylation("Breast - Mammary Tissue", bc_beta_final_BMIQ)
dim(bc_final_BMIQ_breast)
noob_final_BMIQ_breast <- tissue_specific_methylation("Breast - Mammary Tissue", noob_beta_final_BMIQ)
dim(noob_final_BMIQ_breast)
save(raw_final_BMIQ_breast,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/raw_beta_bmiq_2-6-2021/raw_final_BMIQ_breast_2-6-2021.RData")
save(bc_final_BMIQ_breast,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/bc_beta_bmiq_2-6-2021/bc_final_BMIQ_breast_2-6-2021.RData")
save(noob_final_BMIQ_breast,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/noob_beta_bmiq_2-6-2021/noob_final_BMIQ_breast_2-6-2021.RData")

cov_breast <- cov_final_selected[which(cov_final_selected$TSD=="Breast - Mammary Tissue"), ]
cov_breast <- cov_breast[order(cov_breast$CollaboratorSampleID), ]
write.csv(cov_breast, "G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/covariates_tissue_2-6-2021/breast_covariates_2-6-2021.csv", row.names=FALSE)

###Muscle###
raw_final_BMIQ_muscle <- tissue_specific_methylation("Muscle - Skeletal", meth_final_BMIQ)
dim(raw_final_BMIQ_muscle)
bc_final_BMIQ_muscle <- tissue_specific_methylation("Muscle - Skeletal", bc_beta_final_BMIQ)
dim(bc_final_BMIQ_muscle)
noob_final_BMIQ_muscle <- tissue_specific_methylation("Muscle - Skeletal", noob_beta_final_BMIQ)
dim(noob_final_BMIQ_muscle)
save(raw_final_BMIQ_muscle,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/raw_beta_bmiq_2-6-2021/raw_final_BMIQ_muscle_2-6-2021.RData")
save(bc_final_BMIQ_muscle,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019n/bc_beta_bmiq_2-6-2021/bc_final_BMIQ_muscle_2-6-2021.RData")
save(noob_final_BMIQ_muscle,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/noob_beta_bmiq_2-6-2021/noob_final_BMIQ_muscle_2-6-2021.RData")

cov_muscle <- cov_final_selected[which(cov_final_selected$TSD=="Muscle - Skeletal"), ]
cov_muscle <- cov_muscle[order(cov_muscle$CollaboratorSampleID), ]
write.csv(cov_muscle, "G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/covariates_tissue_2-6-2021/muscle_covariates_2-6-2021.csv", row.names=FALSE)

###Kidney###
raw_final_BMIQ_kidney <- tissue_specific_methylation("Kidney - Cortex", meth_final_BMIQ)
dim(raw_final_BMIQ_kidney)
bc_final_BMIQ_kidney <- tissue_specific_methylation("Kidney - Cortex", bc_beta_final_BMIQ)
dim(bc_final_BMIQ_kidney)
noob_final_BMIQ_kidney <- tissue_specific_methylation("Kidney - Cortex", noob_beta_final_BMIQ)
dim(noob_final_BMIQ_kidney)
save(raw_final_BMIQ_kidney,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/raw_beta_bmiq_2-6-2021/raw_final_BMIQ_kidney_2-6-2021.RData")
save(bc_final_BMIQ_kidney,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/bc_beta_bmiq_2-6-2021/bc_final_BMIQ_kidney_2-6-2021.RData")
save(noob_final_BMIQ_kidney,  file="G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/noob_beta_bmiq_2-6-2021/noob_final_BMIQ_kidney_2-6-2021.RData")

cov_kidney <- cov_final_selected[which(cov_final_selected$TSD=="Kidney - Cortex"), ]
cov_kidney <- cov_kidney[order(cov_kidney$CollaboratorSampleID), ]
write.csv(cov_kidney, "G:/Projects/GTEx/GTEX_EPIC_Methylation_September_2019/covariates_tissue_2-6-2021/kidney_covariates_2-6-2021.csv", row.names=FALSE)


