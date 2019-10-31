### R code from vignette source 'Rcade.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: Rcade.Rnw:84-86
###################################################
dir <- file.path(system.file("extdata", package="Rcade"), "STAT1")
dir


###################################################
### code chunk number 3: Rcade.Rnw:92-93
###################################################
library(Rcade)


###################################################
### code chunk number 4: Rcade.Rnw:114-115
###################################################
DE <- read.csv(file.path(dir, "DE.csv"))


###################################################
### code chunk number 5: Rcade.Rnw:120-122
###################################################
DElookup <- list(GeneID="ENSG", logFC="logFC", B="B",
"Genes.Location", "Symbol")


###################################################
### code chunk number 6: Rcade.Rnw:132-133
###################################################
dir(dir, pattern = ".bam")


###################################################
### code chunk number 7: Rcade.Rnw:149-151
###################################################
targets <- read.csv(file.path(dir, "targets.csv"), as.is = TRUE)
targets


###################################################
### code chunk number 8: Rcade.Rnw:165-169
###################################################
anno <- read.csv(file.path(dir, "anno.csv"))

anno <- anno[order(anno$chromosome_name),]
colnames(anno) <- c("ENSG","chr","start","end","str")


###################################################
### code chunk number 9: Rcade.Rnw:174-186 (eval = FALSE)
###################################################
## library(biomaRt)
## 
## anno <- getBM(
## 		attributes= c("ensembl_gene_id", "chromosome_name",
## 			"transcript_start", "transcript_end", "strand"),
## 		mart= useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
## 	)
## 
## ##order, to reduce size of ChIPannoZones object later
## anno <- anno[order(anno$chromosome_name),]
## ##use appropriate column names
## colnames(anno) <- c("ENSG","chr","start","end","str")


###################################################
### code chunk number 10: Rcade.Rnw:191-192
###################################################
ChIPannoZones <- defineBins(anno, zone=c(-1500, 1500), geneID="ENSG")


###################################################
### code chunk number 11: Rcade.Rnw:205-208
###################################################
DE.prior = 0.01
prior.mode = "keepChIP"
prior = c("D|C" = 0.05, "D|notC" = 0.005)


###################################################
### code chunk number 12: Rcade.Rnw:220-222 (eval = FALSE)
###################################################
## library(parallel)
## cl <- makeCluster(4, "SOCK")


###################################################
### code chunk number 13: Rcade.Rnw:227-228
###################################################
cl <- NULL


###################################################
### code chunk number 14: Rcade.Rnw:233-238
###################################################
Rcade <- RcadeAnalysis(DE, ChIPannoZones, annoZoneGeneidName="ENSG",
	ChIPtargets=targets, ChIPfileDir = dir,
	cl=cl, DE.prior=DE.prior, prior.mode=prior.mode, prior=prior,
	DElookup=DElookup)
Rcade


###################################################
### code chunk number 15: Rcade.Rnw:242-243
###################################################
x <- getDE(Rcade)


###################################################
### code chunk number 16: Rcade.Rnw:247-248
###################################################
x <- getChIP(Rcade)


###################################################
### code chunk number 17: Rcade.Rnw:252-253
###################################################
x <- getRcade(Rcade)


###################################################
### code chunk number 18: P1
###################################################
plotPCA(Rcade)


###################################################
### code chunk number 19: P1fig
###################################################
plotPCA(Rcade)


###################################################
### code chunk number 20: P2
###################################################
plotMM(Rcade)


###################################################
### code chunk number 21: P2fig
###################################################
plotMM(Rcade)


###################################################
### code chunk number 22: P3 (eval = FALSE)
###################################################
## library(rgl)
## plotBBB(Rcade)


###################################################
### code chunk number 23: Rcade.Rnw:310-311 (eval = FALSE)
###################################################
## exportRcade(Rcade, directory="RcadeOutput")


###################################################
### code chunk number 24: Rcade.Rnw:316-317 (eval = FALSE)
###################################################
## ?exportRcade


###################################################
### code chunk number 25: Rcade.Rnw:322-323 (eval = FALSE)
###################################################
## exportRcade(Rcade, directory="RcadeOutput", cutoffArg=2000)


###################################################
### code chunk number 26: Rcade.Rnw:327-328 (eval = FALSE)
###################################################
## exportRcade(Rcade, directory="RcadeOutput", cutoffMode="B", cutoffArg=0)


###################################################
### code chunk number 27: Rcade.Rnw:346-347
###################################################
sessionInfo()


