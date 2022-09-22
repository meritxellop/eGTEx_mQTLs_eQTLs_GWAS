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

# DNAm EpiSCORE cells
cells.2=read.csv('data/cell.abundances.csv')

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

EPIC <- read.csv('data/protected/DNAm_GTEx_batch_and_sample_info - MethylationEPIC_GTEX_Batch01_11_SampleSheet.2.csv',header=T)
EPIC$ID <- str_extract(EPIC$Sample_Name,'CO.*')
# Keep sample_plate (1-11), sentrix_id (chip) and sentrix_position (array)
EPIC <- EPIC[,c('ID','Sample_Plate','Sentrix_ID','Sentrix_Position')]

# Path for PEERs
sv.path <- "data/"
sv.files <- list.files(sv.path, full.names = FALSE)
sv.files <-c('BreastMammaryTissue.peers.txt','ColonTransverse.peers.txt','WholeBlood.peers.txt')

# Tissue mappings
tissmap <- read.delim("data/all-tissue-abbreviations.tsv")
tissmap$tissue_id <- gsub('_','',tissmap$tissue_id)

tissmap_2=cbind(c('Breast','Blood','Colon'),c('BreastMammaryTissue','WholeBlood','ColonTransverse'))
colnames(tissmap_2)=c('tiss_ab','tissue_id')
cells.2 <- merge(cells.2,tissmap_2,by.x='Tissue',by.y=1)


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
# For each tissue with surrogate variables, test the cell abundances
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
  print(current.tissue)
  
  # Filter data
  reg.df.a<-subset(cells.2,tissue_id%in%current.tissue)
  for (cell in unique(reg.df.a$Cell)) {
    print(cell)
    reg.df <- subset(reg.df.a,Cell%in%cell)
    rownames(reg.df) <- reg.df$Subject
  
  # Common individuals between metadata and surrogate variables
  common.inds <- intersect(reg.df$Subject, rownames(current.svs))
  reg.df <- reg.df[common.inds, ]
  current.svs <- current.svs[common.inds, ]
  current.svs <- as.matrix(current.svs)
  stopifnot(identical(rownames(current.svs), rownames(reg.df)))
  
  
  # For each cell type, fit a model and get the R2adj
  #for(j in 1:ncol(reg.df)){
    
    tissue <- current.tissue
    varname <- cell
    
    outobj <- tryCatch({
      
      # Fit the model for the current covariate
      X <- reg.df[, 'Fraction']
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
result.df<- merge(result.df,tissmap_2,by.x=1,by.y=2)
colnames(result.df)[2]='Cell'

#############################
# Plotting
#############################


p1 <- ggplot(data=result.df, aes(x=Cell, y=R2adj,fill=tiss_ab)) + geom_bar(stat="identity")+facet_grid(.~tiss_ab,scales='free')+theme_classic()+scale_fill_manual(values=c('#FF00BB','#33CCCC','#EEBB77'))+theme(legend.position = 'none')+labs(tag='a')+ theme(plot.tag = element_text(face='bold'))+ theme(legend.position = "none")+theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white") )
# 1100 x 400

wRPC.o <- read.csv('data/cell.abundances.csv',header=T)

p2 <- ggplot(wRPC.o, aes(x = Cell, y = Fraction, fill = Sex)) + geom_boxplot() + facet_wrap(.~Tissue,scales="free")+ scale_fill_manual(values=c("red", "blue", "grey"))+theme_classic()+ labs(tag = 'b')+theme(plot.tag=element_text(face='bold'))+ theme(legend.position = "none")+theme(legend.key = element_blank(), strip.background = element_rect(colour="white", fill="white") )
# 1100 x 400

plot_grid(nrow = 2,p1,p2)


