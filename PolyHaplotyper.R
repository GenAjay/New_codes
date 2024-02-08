## -----------------------------------------------------------------------------
rm(list=ls()) # clear existing data from memory
library(PolyHaplotyper)
data(PolyHaplotyper_demo) 
# show the demo data:
ls()

## -----------------------------------------------------------------------------
demo_snpdos[1:6, 1:8]

## -----------------------------------------------------------------------------
hblist <- haploblock_df2list(demo_snpdos, mrkcol=1, hbcol=2)
# number of markers in each haploblock:
sapply(hblist, length)

## -----------------------------------------------------------------------------
demo_snpdos <- demo_snpdos[, -2]
demo_snpdos[1:6, 1:8]

## -----------------------------------------------------------------------------
head(demo_ped)

## -----------------------------------------------------------------------------
# create a list of of replicates:
tb <- table(demo_ped$genotype)
replist <- list()
for (dup in names(tb[tb>1])) {
  replist[[dup]] <- demo_ped$sample_nr[demo_ped$genotype == dup]
}

## -----------------------------------------------------------------------------
dim(demo_ped)
for (dup in seq_along(replist)) {
  dupsamp <- as.character(replist[[dup]])[-1]
  demo_ped <- demo_ped[!(demo_ped$sample_nr %in% dupsamp),]
}
dim(demo_ped)

## -----------------------------------------------------------------------------
# merge all the duplicated samples:
dim(demo_snpdos)
demo_snpdos <- mergeReplicates(mrkDosage=demo_snpdos, replist=replist,
                               solveConflicts=TRUE)
dim(demo_snpdos)

## -----------------------------------------------------------------------------
colnames(demo_snpdos) <- 
  demo_ped$genotype[match(colnames(demo_snpdos), demo_ped$sample_nr)]


## -----------------------------------------------------------------------------
parents <- cbind(c(36451, 41234, 9656, 32141),
                 c(39287, 40360, 9541, 39287))
parents
# find the FS individuals by looking in the pedigree for their mother and father:
FS <- list()
for (p in seq_len(nrow(parents))) {
  FS[[p]] <- 
    demo_ped$genotype[!is.na(demo_ped$mother) & demo_ped$mother==parents[p, 1] &
                      !is.na(demo_ped$father) & demo_ped$father==parents[p, 2]]
}
# sizes of the FS families:
sapply(FS, length)

## ----results="hide"-----------------------------------------------------------
results <- inferHaplotypes(mrkDosage=demo_snpdos, ploidy=6,
                           haploblock=hblist,
                           parents=parents, FS=FS)

## -----------------------------------------------------------------------------
names(results[[1]])
# show part of hapdos: the dosages of the haplotypes, for the first haploblock:
results[[1]]$hapdos[, 1:8]
# show the composition of the used haplotypes:
usedhap(results[[1]])
# show the composition of all possible haplotypes:
allhap(results[[1]])

## -----------------------------------------------------------------------------
results[[1]]$imputedGeno[, 1:8]
demo_snpdos[hblist[[1]], colnames(results[[1]]$imputedGeno)[1:8]]

## -----------------------------------------------------------------------------
ovwFS <- overviewByFS(haploblock=hblist, parents=parents, FS=FS,
                      hapresults=results)
# The full table is too wide to show;
# for FS family 1 the results are:
knitr::kable(ovwFS$ovw[, 1: 8], digits=c(0,0,0,0,3,0,0,0))
# the final columns for the non-FS individuals and all individuals:
knitr::kable(ovwFS$ovw[, c(1:2, 27:31)])
#

## -----------------------------------------------------------------------------
pedchk <- pedigreeHapCheck(ped=demo_ped, mrkDosage=demo_snpdos,
                           haploblock=hblist,
                           hapresults=results)
#show part of ped_arr for haploblock 1:
knitr::kable(pedchk$ped_arr[1:8,,1])
#show parents_arr for parents with more than 10 offspring and haploblock 1:
knitr::kable(pedchk$parents_arr[pedchk$parents_arr[,3,1]>10,,1])
#

## -----------------------------------------------------------------------------
cst <- calcStatistics(pedchk=pedchk, ovwFS=ovwFS)
knitr::kable(cst$pedstats)
knitr::kable(cst$FSstats, digits=c(0,0,0,2,2,2))
#

## -----------------------------------------------------------------------------
calcMrkHaptable(ovwFS=ovwFS)

## -----------------------------------------------------------------------------
# show the segregation for FS number 1 (FSnr=1) and 
# haploblock number 1 (hbresults=results[[1]])
showOneFS(FSnr=1, hbresults=results[[1]], mrkDosage=demo_snpdos, 
          FS=FS, parents=parents)

## ----eval=FALSE---------------------------------------------------------------
#  build_ahccompletelist(ploidy=6, maxmrk=5, overwrite=FALSE)

## ----echo=FALSE---------------------------------------------------------------
maxmrk <- 8
maxploidy <- 8
nhapcomb <- matrix(totHapcombCount(ploidy=rep(1:maxploidy, each=maxmrk), 
                                   nmrk=rep(1:maxmrk, maxploidy)),
                   ncol=maxploidy,
                   dimnames=list(nmrk=1:maxmrk, ploidy=1:maxploidy))
nmrkcomb <- matrix((rep(1:maxploidy, each=maxmrk)+1)^rep(1:maxmrk, maxploidy),
                   ncol=maxploidy,
                   dimnames=list(nmrk=1:maxmrk, ploidy=1:maxploidy))
knitr::kable(nhapcomb, caption="haplotype combinations", row.names=TRUE)
knitr::kable(nmrkcomb, caption="combinations of marker dosages", row.names=TRUE)
#

