library(dplyr)
library(plyr)
library(stringr)
library(tibble)
library(tidyr) 
library(ggplot2)
library(ComplexHeatmap)
library(feather)
library(tidyr)
library(readxl)
options(stringsAsFactors = FALSE)

#############################
# Loading data
#############################
# Metadata
individual.meta <- read.delim("data/protected/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt",
                              stringsAsFactors = TRUE)
individual.meta.2 <- read.delim("data/protected/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                              stringsAsFactors = TRUE)
individual.meta <- rbind.fill(individual.meta,individual.meta.2[!individual.meta.2$SUBJID%in%individual.meta$SUBJID,])

sample.meta <- read.delim("data/protected/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt",
                              stringsAsFactors = TRUE,row.names=1)

# Metadata description
individual.meta.desc <- read_excel("data/protected/GTEx_Analysis_2016-01-15_v7_Annotations_SubjectPhenotypesDD.xlsx") %>%
    select(VARNAME, VARDESC)

# xCell scores, not normalized
cells=as.data.frame(t(read.table('data/GTEx_Analysis_v8_xCell_scores_7_celltypes.csv',sep=',',header=T,check.names=F,row.names=1)))
# Select cells analyzed in Kim-Hellmuth et al. 2020
cells <- cells[,c('Adipocytes','Epithelial_cells','Neutrophils','Myocytes')]
cells$Tissue <- sample.meta[rownames(cells),'SMTSD']
cells$Subject <- str_extract(rownames(cells),'GTEX-\\w+')
cells$Myocytes[!cells$Tissue %in% 'Muscle - Skeletal'] <- NA 
cells$Neutrophils[!cells$Tissue %in% 'Whole Blood'] <- NA 
cells$Adipocytes[!cells$Tissue %in% 'Breast - Mammary Tissue'] <- NA 
cells$Epithelial_cells[cells$Tissue %in% c('Whole Blood','Ovary','Muscle - Skeletal', 'Testis')] <- NA 

vardesc <- individual.meta.desc
vardesc=as.data.frame(cbind(vardesc$VARNAME,vardesc$VARDESC))
vardesc <- rbind(vardesc,cbind(c('DTHSEASON','Sample_Plate','Sentrix_ID','Sentrix_Position','Autolysis.Score','CigSmoke_NFC'),c('Death Season','Plate','Chip','Array','Autolysis Score','Smoking Status')))
#vardesc <- rbind(vardesc,cbind(colnames(cells)[1:4],colnames(cells)[1:4]))
colnames(vardesc)=c('VARNAME','VARDESC')

# Fix metadata error
vardesc[vardesc$VARNAME == "HGHT", "VARDESC"] <- "Height"

# Manual fixes
vardesc[vardesc$VARNAME == "MHT2D", "VARDESC"] <- "Diabetes mellitus type II"
vardesc$VARDESC[vardesc$VARNAME%in%'MHCOPD'] <- gsub('\\(.*','',vardesc$VARDESC[vardesc$VARNAME%in%'MHCOPD'])
vardesc$VARDESC[vardesc$VARNAME%in%'MHHRTDIS'] <- gsub('\\(.*','',vardesc$VARDESC[vardesc$VARNAME%in%'MHHRTDIS'])
vardesc$VARDESC[vardesc$VARNAME%in%'SEX'] <- 'Biological sex (male or female)'

smoking <- read.csv('data/protected/gtex_smoking_4-22-2021.csv',header=T)
smoking$CigSmoke_NFC <- as.numeric(gsub('Former',1,gsub('Never',0,gsub('Current',2,smoking$CigSmoke_NFC))))
smoking <- smoking[,c('CollaboratorParticipantID','TSD','CigSmoke_NFC')]

eGTEx <- read.csv('data/protected/DNAm_GTEx_batch_and_sample_info - eGTExDNA_Pierce_Jan\'18.csv',header=T)
eGTEx$ID <- paste0(eGTEx$Container,'_',eGTEx$Position)
eGTEx <- eGTEx[,c('ID','Collaborator.Participant.ID','Tissue.Site.Detail','Autolysis.Score')]
eGTEx$Autolysis.Score <- as.numeric(gsub('Severe',3,gsub('Moderate',2,gsub('None',0,gsub('Mild',1,eGTEx$Autolysis.Score)))))

#EPIC <- read.csv('data/DNAm_GTEx_batch_and_sample_info - MethylationEPIC_GTEX_Batch01_11_SampleSheet.csv',header=T,skip=7)
EPIC <- read.csv('data/protected/DNAm_GTEx_batch_and_sample_info - MethylationEPIC_GTEX_Batch01_11_SampleSheet.2.csv',header=T)
EPIC$ID <- str_extract(EPIC$Sample_Name,'CO.*')
# Keep sample_plate (1-11), sentrix_id (chip) and sentrix_position (array)
EPIC <- EPIC[,c('ID','Sample_Plate','Sentrix_ID','Sentrix_Position')]

# Path for PEERs
sv.path <- "data/"
sv.files <- list.files(sv.path, full.names = FALSE)
sv.files <-c('BreastMammaryTissue.peers.txt','ColonTransverse.peers.txt','KidneyCortex.peers.txt','Lung.peers.txt','MuscleSkeletal.peers.txt','Ovary.peers.txt','Prostate.peers.txt','Testis.peers.txt','WholeBlood.peers.txt')

# Tissue mappings
tissmap <- read.delim("data/all-tissue-abbreviations.tsv")
tissmap$tissue_id <- gsub('_','',tissmap$tissue_id)

#############################
# Preprocessing
#############################
# Prepare all metadata (individual + sample) for samples with expression

samp.withdata <- merge(eGTEx,smoking,by.x=c('Collaborator.Participant.ID','Tissue.Site.Detail'),by.y=c('CollaboratorParticipantID','TSD'))
samp.withdata <- merge(samp.withdata,individual.meta,by.x='Collaborator.Participant.ID',by.y='SUBJID');
samp.withdata <- merge(samp.withdata,EPIC,by.x='ID',by.y='ID');
samp.withdata <- merge(samp.withdata,tissmap[,1:2],by.x='Tissue.Site.Detail',by.y='tissue_site_detail')
#samp.withdata <- merge(samp.withdata,cells,by.x=c('Collaborator.Participant.ID','Tissue.Site.Detail'),by.y=c('Subject','Tissue'),all.x=T)

# Medical history variables are factors
mhvars <- colnames(samp.withdata)[grep("MH.*", colnames(samp.withdata))]
samp.withdata[, mhvars] <- lapply(samp.withdata[, mhvars], as.factor)

# String columns to exclude
toexclude <- c("MHGENCMT", "DTHLUCOD", "SMPTHNTS", "SMTS", "SMTSC", "DTHTIME", "MHTTCMT", "MHBLDDNDR")
toexclude <- toexclude[toexclude%in%colnames(samp.withdata)]
samp.withdata <- samp.withdata %>%
    select(-toexclude)

# Other variables to turn into factors
tofactors <- "Sample_Plate Sentrix_ID Sentrix_Position LBCMVTAB LBEBVGAB LBEBVMAB LBHBCABM LBHBCABT LBHBSAB LBHBSAG LBHCV1NT LBHBHCVAB LBHIV1NT LBHIVAB LBHIVO LBPRRVDRL LBRPR SEX RACE ETHNCTY"
tofactors <- strsplit(tofactors, split = " ")[[1]]
tofactors <- tofactors[tofactors%in%colnames(samp.withdata)]
samp.withdata[, tofactors] <- lapply(samp.withdata[, tofactors], as.factor)

# Remove variables with many levels (provided by Diego) or other
# variables that do not make sense to test
many.level.ones <- c("SMPTHNTS", "TRCHSTIN", "DTHCOD", "DTHFUCOD", "SMNABTCH",
                     "SMNABTCHD", "SMGEBTCH", "SMGEBTCHD", "DTHPRNINT", "TRCRTMPL", "TRCRTMPU",
                     "TRISCH")
many.level.ones <- many.level.ones[many.level.ones%in%colnames(samp.withdata)]
samp.withdata <- samp.withdata[, !colnames(samp.withdata) %in% many.level.ones]


#############################
# Analysis
#############################
# For each tissue with surrogate variables, test the covariates
result.df <- data.frame(Tissue = character(),
                        Variable = character(),
                        R2adj = numeric(),
                        Warn = character(),
                        stringsAsFactors = FALSE)
k <- 1

for(i in seq_along(sv.files)){

    # Read PEERs, keep top 3
    current.sv.file <- sv.files[i]
    current.svs <- read.table(paste0(sv.path, current.sv.file))[,1:3]

    # Get name for current tissue
    current.tissue <- gsub("(.*)(.peers.txt)", "\\1", current.sv.file)

    # Filter data
    reg.df <- samp.withdata %>%
        filter(tissue_id == current.tissue) %>%
        column_to_rownames("Collaborator.Participant.ID") %>%
        select(-ID, -Tissue.Site.Detail, -tissue_id)
    
    # Common individuals between metadata and surrogate variables
    common.inds <- intersect(rownames(reg.df), rownames(current.svs))
    reg.df <- reg.df[common.inds, ]
    current.svs <- current.svs[common.inds, ]
    current.svs <- as.matrix(current.svs)
    stopifnot(identical(rownames(current.svs), rownames(reg.df)))

    # Sanity checks
    # Keep variables with more than 80% of cases, but do not apply such filter for cell abundances
    #reg.df.cells <- reg.df[,colnames(reg.df)%in%colnames(cells)]
    #reg.df <- reg.df[,!colnames(reg.df)%in%colnames(cells)]
    reg.df <- reg.df[, apply(reg.df, 2, function(x){sum(!is.na(x))})/nrow(reg.df) >= 0.1]
    #reg.df <- cbind(reg.df, reg.df.cells)

    # If it's a medical history variable, require at least 20% of cases
    mhcheck <- reg.df[, grep("MH.*", colnames(reg.df))]
    mhcheck <- mhcheck[, !colnames(mhcheck) %in% c("MHSRC")]
    mhcheck <- apply(mhcheck, 2, function(x){sum(as.numeric(as.character(x)), na.rm=T)})/nrow(mhcheck)
    mhcheck <- names(mhcheck[mhcheck <= 0.2])
    reg.df <- reg.df[, !colnames(reg.df) %in% mhcheck]
    
    # For each covariate, fit a model and get the R2adj
    for(j in 1:ncol(reg.df)){

        tissue <- current.tissue
        varname <- colnames(reg.df)[j]
        
        outobj <- tryCatch({

            # Fit the model for the current covariate
            X <- reg.df[, j]
            Y <- current.svs
            #Y <- current.svs[,c(1,1)] # first PEER only
                
            # Rows with nas, remove
            naidx <- is.na(X)
            Y <- Y[!naidx, ]
            X <- X[!naidx]
            
            fitted.model <- lm(Y ~ X, na.action = "na.omit")

            ## # Calculating multivariate R2
            ## # http://users.stat.umn.edu/~helwig/notes/mvlr-Notes.pdf
            ## ybar <- colMeans(Y)
            ## n <- nrow(Y)
            ## m <- ncol(Y)
            ## Ybar <- matrix(ybar, n, m, byrow=TRUE)
            ## SSCP.T <- crossprod(Y - Ybar)
            ## SSCP.R <- crossprod(fitted.model$fitted.values - Ybar)
            ## SSCP.E <- crossprod(Y - fitted.model$fitted.values)

            ## # Degrees of freedom
            ## m <- ncol(Y)
            ## n <- nrow(Y)
            ## p <- 1
            
            ## dfT <- m*(n-1)
            ## dfR <- m*p
            ## dfE <- m*(n-p-1)
            
            ## R2 <- sum(diag(SSCP.R))/sum(diag(SSCP.T))
            ## R2adj <- 1 - (sum(diag(SSCP.E))/dfE)/(sum(diag(SSCP.T))/dfT)

            ### Get the mean adjusted Rsquared
	    ### PEER1 should weight more than PEERn. Weight average R2adj across PEERs by PEER reciprocal rank: http://www.gitta.info/Suitability/en/html/Normalisatio_learningObject1.html
	    
	    m <- ncol(Y) # Number of PEERs
            weights_per_peer <- unlist(lapply(seq(1,m),FUN = function(x) {1/x}))
            R2adj <- lapply(X=summary(fitted.model), FUN=function(x){ return(x$adj.r.squared) })
            R2adj <- weighted.mean(unlist(R2adj),weights_per_peer)
	    if (sign(R2adj) < 0) { R2adj <-0 }
            
            #R2adj <- mean(unlist(R2adj))
            
            # If no failures, return calculated values
            list(R2adj = R2adj, Warn = NA)
            
        }, warning = function(w){

            msg <- conditionMessage(w)
            
            return(list(R2adj = NA, Warn = msg))
            
        }, error = function(e){

            msg <- conditionMessage(e)

            return(list(R2adj = NA, Warn = msg))
        })

        # Saving data
        result.df[k, "Tissue"] <- tissue
        result.df[k, "Variable"] <- varname
        result.df[k, "R2adj"] <- outobj$R2adj
        result.df[k, "Warn"] <- outobj$Warn

        k <- k +1
    }
}

result.df <- result.df %>%
    mutate(R2adj = ifelse(is.nan(R2adj), NA, R2adj))

#############################
# Plotting
#############################
plotdf <- result.df %>%
    select(-Warn) %>%
    pivot_wider(names_from = c("Variable"), values_from = "R2adj") %>%
    as.data.frame %>%
    column_to_rownames("Tissue")
    
# Keep variables with an R2 larger than 0.02 in at least one tissue
vars.keep <- apply(plotdf, 2, function(x){sum(x >= 0.02, na.rm=T) >= 1})
vars.keep <- vars.keep[vars.keep]
vars.keep <- names(vars.keep)
plotdf <- plotdf[, vars.keep]

# If it's NA In more than 5 tissues, remove
#plotdf.c <- plotdf[,c('Epithelial_cells','Neutrophils')]
#plotdf <- plotdf[,!colnames(plotdf)%in%c('Epithelial_cells','Neutrophils')]
#plotdf <- cbind(plotdf[ , apply(plotdf, 2, function(x){sum(is.na(x))}) <= 5],plotdf.c)
plotdf <- plotdf[ , apply(plotdf, 2, function(x){sum(is.na(x))}) <= 5]

# Replace variable names by descriptions
plotcols <- data.frame(VARNAME = colnames(plotdf),
                       stringsAsFactors = FALSE) %>%
    left_join(x = .,
              y = vardesc,
              by = "VARNAME")
plotcols$VARDESC <- gsub("(.*)(:.*$)", "\\1", plotcols$VARDESC)

colnames(plotdf) <- plotcols$VARDESC

# Replace long tissue names by abbreviations
#rownames(plotdf) <- tissmap[match(rownames(plotdf), tissmap$tissue_id), "tissue_abbrv"]

# Reorder the matrix according to hclust
roword <- hclust(dist(plotdf))$order
plotdf2 <- plotdf
# Dirty trick: impute WholeBlood-variable NAs, to keep Neutrophil column and avoid NAs when calculing dist
#plotdf2['WholeBlood',is.na(plotdf2['WholeBlood',])] <- t(as.data.frame(colMeans(plotdf[,is.na(plotdf['WholeBlood',])],na.rm=T)))
colord <- hclust(dist(t(plotdf2)))$order

# For the column order, we bring the one corresponding to "Plate" to the very end
# so we can split the heatmap and highlight it
valbottom <- colord[which(colnames(plotdf[, colord]) %in% c("Plate","Array","Chip"))]
colord <- c(colord[-which(colord %in% valbottom)], valbottom)
plotdf <- plotdf[roword, colord]
plotdf <- t(plotdf)

# MANUAL EDIT!!!!
tissue_abb=c('Blood','Breast','Muscle','Colon','Lung','Ovary','Testis','Kidney','Prostate')
names(tissue_abb)=colnames(plotdf)

# Build circle annotation
ann.df <- data.frame(tissue_id = colnames(plotdf),
                     stringsAsFactors = FALSE) %>%
    left_join(x = .,
              y = tissmap %>% select(tissue_id, color_hex, tissue_site_detail),
              by = "tissue_id") %>%
    mutate(color_hex = paste0("#", color_hex))    

## tissue.colors <- data.frame(tissue_abbrv = colnames(part2)) %>%
##     left_join(x = .,
##               y = tissmap %>% select(tissue_abbrv, color_hex),
##               by = "tissue_abbrv") %>%
##     mutate(color_hex = paste0("#", color_hex))

## tissue.colors <- setNames(tissue.colors$color_hex, tissue.colors$tissue_abbrv)
tissue.colors <- setNames(ann.df$color_hex, ann.df$tissue_abbrv)
colpal <- setNames(rep("#000000", nrow(ann.df)), ann.df$Abbrev)

#colnames(plotdf) <- tissue_abb[colnames(plotdf)]
p1 <- Heatmap(plotdf,
        col = rev(c("#08509A", "#5198e0","#F7FBFF")),
        na_col = "grey",
	row_names_max_width = unit(20, "cm"),
#        bottom_annotation = ha,
        heatmap_legend_param = list(title = expression(R[adj]**2),
                                    at = c(0, 0.1, 0.20, 0.35)),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_split = c(rep(1, nrow(plotdf)-3), 2, 2, 2), # Corresponding to Plate, Array and Chip variables
        row_title = c("", ""),
        border = TRUE)

# Tissue colors
colnames(plotdf) <- c('Colon','Blood','Testis','Kidney','Prostate','Breast','Muscle','Lung','Ovary')
ha <- columnAnnotation(Tissue = anno_simple(x=rep(NA, length(tissue.colors)),
                                           na_col = "white",
                                           pch = 16,
                                           pt_gp = gpar(col = tissue.colors),
                                           pt_size = unit(5, "mm"),
                                           ),
                      show_annotation_name=FALSE)

#pdf("vw_sv-vs-covariates.no_cells.pdf", width = 10, height = 11)
draw(ha %v% p1)
#dev.off()
