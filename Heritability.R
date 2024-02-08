
rm(list=ls())
library(sommer)
library(SpATS)
library(lme4)
library(car)


### Heritability for each model
setwd("")
# chech dimensions pheno = geno
pheno.v0 <- read.csv("pheno.csv")
pheno.v0$Whole.milling.percentage.LSU_01.0000100
# Genotype -------------------------------------------------------------------
# ~ GBLUP: GRM
geno1 <- readRDS("geno.rds")
geno2 <- readRDS("geno2.rds") # in case you have two files
geno <- rbind(geno1, geno2)
geno.v <- geno[rownames(geno) %in% pheno.v0$germplasmName, ]

# efectos aditivos, de dominancia y espistasis
G_mat <- A.mat(X=as.matrix(geno.v), min.MAF = .05)
#add the other matrices with sommer functions

#### Model 1: additive effects - Calculating h2 ####------------------------------------------

pheno.v1 <- pheno.v0[pheno.v0$germplasmName %in% rownames(G_mat), ]
#pheno.v1$germplasmName2 <- pheno.v1$germplasmName # epistasis or dominance accounted with a cloned id factor (if you need)

pheno.v1$rowf <- as.factor(pheno.v1$rowNumber)
pheno.v1$colf <- as.factor(pheno.v1$colNumber)
pheno.v1$repf <- as.factor(pheno.v1$replicate)

mm.a <- mmer(fixed = yield ~ 1, #lets do for yield but here you put your trait
             random = ~vsr(germplasmName,Gu=G_mat),
             rcov = ~ vsr(units),
             data = pheno.v1
)

#summary(mm.a)$varcomp
# heritability
#vpredict(mm.a, h2 ~ (V1) / ( V1+V2)) # broad sense
vpredict(mm.a, h2 ~ (V1) / ( V1+V2)) # narrow sense

#Spatial correction
mm.sp.a <- mmer(fixed = yield ~ 1, 
                random = ~vsr(germplasmName,Gu=G_mat) +
                  vsr(rowf) + vsr(colf) + spl2Da(rowNumber,colNumber),
                rcov = ~ vsr(units),
                data = pheno.v1
)

#summary(mm.sp.a)$varcomp

# heritability
#vpredict(mm.sp.a, h2 ~ (V1) / ( V1+V5)) # broad sense
vpredict(mm.sp.a, h2 ~ (V1) / ( V1+V5)) # narrow sense

############################################################################
rm(list=ls())
library(sommer)
library(SpATS)
library(lme4)
library(car)

# Reading data --------------------------------------------------------------
setwd("")

# Phenotype -----------------------------------------------------------------
pheno.v0 <- read.csv("pheno", header = T)
#pheno.v0$Grain.length.average.LSU_01.0000102

# checking outliers from the pheno ------------------------------------------

# outlier detection and elimination
fit <- lm(yield ~ germplasmName, data = pheno.v0)
(outlier <- names(outlierTest(fit)$p))
pheno.v0[outlier, "yield"] <- NA



# Adjust phenotype by the spatial trends
pheno.v1 <- pheno.v0

pheno.v1$rowf <- as.factor(pheno.v1$rowNumber)
pheno.v1$colf <- as.factor(pheno.v1$colNumber)
pheno.v1$repf <- as.factor(pheno.v1$replicate)

nrow <- max(pheno.v1$rowNumber)
ncol <- max(pheno.v1$colNumber)
nseg.row <- nrow
nseg.col <- ncol

#this is a loop to do for all traits at once
trait_columns <- c("yield", "dth", "milling", 
                   "chalk")#,
#"grain_length",


# Initialize a list to store adjusted results for each trait
adjusted_results <- list()

# Loop through each trait
for (trait in trait_columns) {
  # weights for ID's - adjust residual for further analysis 
  
  fitR <- SpATS(response = trait, 
                fixed = ~ 1, 
                random = ~ colf + rowf,# + repf, 
                spatial = ~ PSANOVA(colNumber, rowNumber, nseg = c(nseg.col, nseg.row)), 
                genotype = "germplasmName", 
                genotype.as.random = TRUE, 
                data = pheno.v1)
  
  
  vcov.mme <- fitR$vcov$C11_inv 
  w <- diag(vcov.mme)
  # extracting the blues
  blups <- SpATS::predict.SpATS(fitR,which = "germplasmName" )
  adjusted <- cbind(blups, w)
  
  adjusted_results[[trait]] <- adjusted
}

# Combine results from the list into a single data frame
combined_results <- do.call(rbind, lapply(names(adjusted_results), function(trait) {
  data.frame(germplasmName = rownames(adjusted_results[[trait]]), trait = trait, adjusted_value = adjusted_results[[trait]])
}))

write.csv(combined_results, "blues.csv")

############################################################################Multiyear adjustments########
rm(list=ls())
library(sommer)
library(SpATS)

##### set working directory###

setwd("")

##### Loading the data#####


blup1 <- read.csv("trial1_blues.csv", header=T)
blup2 <- read.csv("trial2_blues.csv", header=T)
blup<- rbind(blup1, blup2)

##### MET Model 

blup$trial <- as.factor(blup$trial)
blup$germplasmName <- as.factor(blup$germplasmName)

str(blup) #make sure your traits are numerics


# Calculate the BLUPS
met.<- mmer  (fixed = milling ~ 1 + trial,  
              random = ~ germplasmName,
              rcov = ~ vsr(units),
              data = blup)

# with weights
met.<- mmer  (fixed = milling ~ 1 + trial,  
              random = ~ germplasmName,
              rcov = ~ vsr(units),
              data = blup,
              weights = milling_w )


# One trial, to have just blups (one date of planting)

met.<- mmer  (fixed = chalk ~ 1,  
              random = ~ germplasmName,
              rcov = ~ vsr(units),
              data = blup,
              weights = chalk_w )

# extract the blups form any of these models
blups.met <- predict.mmer(met., classify=c("germplasmName"), D = "germplasmName" )
blups <- blups.met$pvals

write.csv(blups, "blups.csv" )

##### Loading the data for making a loop for all the traits#####


blup1 <- read.csv("trial1.csv", header=T)
blup2 <- read.csv("trial2.csv", header=T)

selected_columns <- c("yield", "dth", "grain_length", "chalk", "milling", "trial", "germplasmName")  
data1_subset <- blup1[, selected_columns] # when both data sets don't have same columns
data2_subset <- blup2[, selected_columns]
blup <- rbind(data1_subset, data2_subset)

blup<- rbind(blup1, blup2)
str(blup)

##### MET Model 

blup$trial <- as.factor(blup$trial)
blup$germplasmName <- as.factor(blup$germplasmName)
#blup$milling <- as.numeric(blup$milling)

### loop for all the traits, with no weights

trait_columns <- c("yield", "dth", "grain_length", "chalk", "milling")


# Initialize a list to store adjusted results for each trait
adjusted_results <- list()

# Loop through each trait
for (trait in trait_columns) {
  
  met_model <- mmer(fixed = as.formula(paste(trait, "~ 1 + trial")),  
                    random = ~ germplasmName, # Random effect of germplasmName
                    rcov =  ~ vsr(units),    
                    data = blup)
  
  # Calculate multi-environmental BLUPs
  blups_met <- predict.mmer(met_model, classify = "germplasmName", D = "germplasmName")
  blups <- blups_met$pvals
  
  adjusted_results[[trait]] <- blups
}

# Combine results from the list into a single data frame
combined_results <- do.call(rbind, lapply(names(adjusted_results), function(trait) {
  data.frame(germplasmName = rownames(adjusted_results[[trait]]), trait = trait, adjusted_value = adjusted_results[[trait]])
}))

write.csv(combined_results, "blups.csv")


#----FIN----#

##########################################################Prediction######################################

rm(list=ls())
library(sommer)
library(ggplot2)

##### Working Directory and Data Loading ####

setwd("")

#### Training set
pheno.tr <- read.csv("Training.csv",header = T)

#### Validation Set
pheno.v0 <- read.csv("prediction.csv",header = T)
pheno.tr2 <- pheno.tr[!pheno.tr$germplasmName %in% pheno.v0$germplasmName,]

### All Phenotype
pheno.v0$chalk=NA
pheno.all <- rbind(pheno.tr2,pheno.v0)
str(pheno.v0)

#### Genotype and GRM
geno <- readRDS("geno.last.rds")
geno.v <- geno[rownames(geno) %in% pheno.v0$germplasmName, ]
geno.tr <- geno[rownames(geno) %in% pheno.tr2$germplasmName, ]
geno.all <- rbind(geno.v,geno.tr)
G_mat <- A.mat(X=as.matrix(geno.all), min.MAF = .05)

# Model using only additive effects
mm.gp <- mmer(fixed = chalk ~ 1, 
              random = ~vsr(germplasmName,Gu=G_mat),
              rcov = ~ vsr(units),
              data = droplevels(pheno.all[pheno.all$germplasmName %in% rownames(geno.all),])
)

pheno.v0 <- read.csv("prediction.csv",header = T)
gebv.blup <- merge(pheno.v0, 
                   data.frame(germplasmName=names(mm.gp$U$`u:germplasmName`$chalk),
                              gebv=as.numeric(mm.gp$U$`u:germplasmName`$chalk)),
                   by = "germplasmName")

cor(gebv.blup$gebv, gebv.blup$chalk)

# ---- Fin----

###############Crossvalidation k fold##############
rm(list=ls())
library(sommer)
library(tidyverse)
library(ggplot2)
library(paletteer)
library(hrbrthemes)

setwd("~/LSU Maria Montiel/Projects/project_mp/project_mp_01paper/project_mp_epistasis/Cross validation/data/")
my_palette <- c("#EE3A8C", "#4682B4", "#2F4F4F","#FFA500")


######------------------Cross Validation K Fold-----------------------######


pheno <- read.csv("blues_corrected_spatial.csv")
geno <- readRDS("geno.last.rds")
geno.v <- geno[rownames(geno) %in% pheno$germplasmName, ]
G_mat <- A.mat(X=as.matrix(geno.v), min.MAF = .05)
G_e_mat <- E.mat(X=as.matrix(geno.v), min.MAF = .05)
pheno$germplasmName2 <- pheno$germplasmName



ID <- colnames(pheno)[1]
res1 <- data.frame()
cyc=25
k=5

for (r in 1:25){
  df. <- pheno
  df. <- df.[sample(1:nrow(df.),nrow(df.)),]
  G_mat_CV <- G_mat
  G_mat.v <- G_mat
  pred.loo <- data.frame()
  for (tr in c("yield", "milling", "dth", "chalk"#, "grain_length"
  )){
    acc.k=NULL
    flds <- caret::createFolds(seq(1:nrow(df.)),
                               k = k, list = TRUE, returnTrain = FALSE)
    pred.i <- data.frame()
    for(fd in 1:k){
      df <- df.
      pred.k <- data.frame()
      geno.k <- df[flds[[fd]],"germplasmName"]
      df[df[,ID]%in%geno.k,tr] <- NA    #NA masking of the lines in the k-fold to be predicted
      df[,ID] <- as.factor(df[,ID])
      frm.cv <- paste(tr,"~ 1  ",collapse = " ")
      fix_frm.cv <- formula(paste(frm.cv, collapse = " "))
      mm.gb <-mmer(fixed = fix_frm.cv,
                   random =  ~  vsr(germplasmName ,Gu=G_mat.v) ,
                   rcov = ~vsr(units),
                   #tolParInv = 100000,
                   data = droplevels(df[df$germplasmName %in% rownames(G_mat.v),])
      )
      pred <- data.frame(germplasmName=names(mm.gb$U[[1]][[tr ]]),
                         gebv=as.numeric(mm.gb$U[[1]][[tr ]]))
      
      rownames(pred) <- pred[,1]
      pred.k <- rbind(pred.k,pred[geno.k ,])
      pred <- merge(pred.k,df.[,c(ID,tr)],by=ID)
      acc.k <- c(acc.k,cor(pred[,2],pred[,3],use = "complete.obs"))
      pred.i <- rbind(pred.i,pred.k)
    }
    pred.i <- merge(pred.i,df.[,c(ID,tr)],by=ID)
    
    res1 <- rbind(res1,data.frame(#pop_id=pop,
      trait=tr,n=nrow(df.),
      acc.k=mean(acc.k,na.rm = T),
      acc.i=cor(pred.i[,2],pred.i[,3],use = "complete.obs"))
    )
    
  }
}

write.csv(res1, "CV.csv")

## fin
###################
rm(list=ls())
library(sommer)
library(tidyverse)
library(ggplot2)
library(paletteer)
library(hrbrthemes)

## Figure for cross validation

setwd("")
my_palette <- c("#EE3A8C", "#2F4F4F","#4682B4","#FFA500")

res <- read.csv("CV.csv")

ggplot(res, aes(x=Material, y=Accuracy, fill=Model))+
  #geom_jitter(color="black", size=0.6, alpha=0.9)+
  geom_boxplot(color="black",position=position_dodge(preserve = "single"))+
  scale_fill_manual(values=my_palette)+
  facet_wrap(~trait)+
  theme_bw(14)+
  ggtitle("Cross Validation Predictive Ability of your pops") +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size=16)) + 
  labs(y = "Predictive ability", x = "Materials") +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 1, by = 0.2))



# with the glitter
ggplot(res, aes(x=trait, y=Accuracy, fill=Material))+
  geom_jitter(color="black", size=0.6, alpha=0.9)+
  geom_boxplot(width=0.40, alpha= 0.7)+ 
  scale_fill_brewer(palette = "Dark2") +
  theme_ipsum() +
  theme (axis.title.x = element_text(size = 12),
         axis.text.x = element_text(size = 10),
         axis.title.y = element_text(size = 10), 
         legend.position = "none")

#######################################################################################







