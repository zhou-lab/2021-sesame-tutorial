## Installing SeSAMe:
## Scripts used in this video:
install.packages(“BiocManager”)
BiocManager::install(version=“devel”)
BiocManager::install(“sesame”)
library(sesame)
sesameDataCacheAll()
packageVersion(“sesame”)

## For more information, visit https://www.bioconductor.org/packages/release/bioc/html/sesame.html

## Pre-processing
## Scripts used in this video:
dest_dir = tempdir()
dest_dir
setwd(dest_dir)
## Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2178224
download.file(“https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2178nnn/GSM2178224/suppl/GSM2178224_184AA3_Grn.idat.gz”, “GSM2178224_184AA3_Grn.idat.gz”)
download.file(“https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2178nnn/GSM2178224/suppl/GSM2178224_184AA3_Red.idat.gz”, “GSM2178224_184AA3_Red.idat.gz”)
download.file(“https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2178nnn/GSM2178225/suppl/GSM2178225_184AA2_Grn.idat.gz”, “GSM2178225_184AA2_Grn.idat.gz”)
download.file(“https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2178nnn/GSM2178225/suppl/GSM2178225_184AA2_Red.idat.gz”, “GSM2178225_184AA2_Red.idat.gz”
list.files(pattern-“*.idat”)
library(sesame)
s = readIDATpair(“GSM2178224_184AA3”)
s
mft = sesameDataGet(“EPIC.address”)$ordering
head(mft)
readIDATpair(“GSM2178225_184AA2”, manifest=mft)
library(parallel)
mclapply(searchIDATprefixes(“.”), readIDATpair, mc.cores=2)
sesameQC(s)
qualityRank(s)$rank_probe_success_rate
betas = getBetas(s)
head(betas, 20)
s
head(s$mask)
sum(s$mask)
s0=resetMask(s)
sum(s0$mask)
s1= qualityMask(s0)
s2=pOOBAH(s1)
sum(s2$mask)
pval = pOOBAH(s, return.pval=TRUE)
s3 = addMask(s1, pval>0.05)
sum(s3$mask) == sum(s2$mask)
sesamePlotIntensVsBeta(s)
s4 = noob(s)
sesamePlotIntensVsBetas(s4)
sesamePlotRedGrnQQ(s4)
s5 = dyeBiasNL(s4)
sesamePlotRedGrnQQ(s5)
sesamePlotIntensVsBetas(s5)
betas=do.call(cbind,mclapply(searchIDATprefixes(“.”), function(px) getBetas(dyeBiasNL(noob(pOOBAH(readIDATpair(px))))), mc.cores=2))
head(betas)
betas2 = openSesame(“.”)
all(betas == betas2, na.rm=T)

## For more information, visit https://www.bioconductor.org/packages/release/bioc/html/sesame.html

## Modeling Differential Methylation
## Scripts used in this video:
library(sesame)
library(SummarizedExperiment)
library(tidyverse)
packageVersion(“sesame”)
se = sesameDataGet(“MM285.20Kx467.SE”)
meta=dplyr::select(as_tibble(colData(se)), IDAT, Sex, Age, Strain, Tissue)
meta
meta$Sex = relevel(factor(meta$Sex), ref=“Female”)
meta$Strain= relevel(factor(meta$Strain), ref=“C57BL/6J”)
meta$Tissue=relevel(factor(meta$Tissue), ref=“Colon”)
str(meta)
betas = assay(se)
ok1 = checkLevels(beats,meta$Sex)
sum(ok1)
betas[which(!ok1), [1], ]
ok2 = checkLevels(betas,meta$Strain)
ok3 = checkLevels(betas,meta$Tissue)
betas = betas[ok1&ok2&ok3, ]
dim(betas)
smry = DML(betas, ~Sex+Age+Strain+Tissue, meta=meta, mc.cores=4)
smry
smry[[1]]
res = summaryExtractTest(smry)
dim(res)
colnames(res)
res %>% arrange(Est_Age) %>% select(Est_Age, Pval_Age) %>% tail
res %>% arrange(Est_Age) %>% select(Est_Age, Pval_Age) %>% head
ggplot(tibble(betaValue = betas[“cg29499259_BC21”, ], age = meta$Age), aes(age,betaValue)) + geom_point() + geom_smooth(method=”lm”)
res %>% arrange(Est_SexMale) %>% select(Est_SexMale, Pval_SexMale) %>% head()
with(res,plot(Est_SexMale, -log10(Pval_SexMale), xlab=“Delta Beta”, ylab=“-log10(P-value)”))
 res %>% filter(Est_SexMale > 0.1, Pval_SexMale < 0.01) %>% rownames_to_column(“Probe_ID”) %>% attachManifest() %>% with(table(seqnames()))
 res %>% filter(Est_SexMale > 0.1, Pval_SexMale < 0.01) %>% rownames_to_column(“Probe_ID”) %>% attachManifest() %>% with(table(seqnames))
res %>% filter(Eff_Tissue > 0.1, FPval_Tissue < 0.01) %>% arrange(-Eff_Tissue) %>% select(Eff_Tissue, FPval_Tissue) %>% head()
withref = res %>% select(starts_with(“Est_Tissue”)) %>% mutate(Est_TissueColon=0)
apply(withref – apply(withref,1,median),2,function(x) sum(x<-0.3))
barplot(apply(withref–apply(withref,1,median),2,function(x) sum(x< -0.3)), las=2)
merged = DMR(betas, smry, “TissueStomach”)
head(merged)
merged = DMR(betas,smry, “SexMale”)
head(merged)

## For more information, visit https://www.bioconductor.org/packages/release/bioc/html/sesame.html

## Inferring MetaData
## Scripts used in this video:
library(sesame)
library(tidyverse)
sh=sesameDataGet(“EPIC.1.SigDF”)
sm=sesameDataGet(“MM285.1.SigDF”)
sh
sm
inferSex(sh)
pOOBAH(sh) %>% attachManifest %>% dplyr::filter(seqnames=“chrY”) %>% with(sum(mask) / length(mask))
inferSexKaryotypes(sh)
inferEthnicity(sh)
inferStrain(getBetas(sm, mask=FALSE))
sort(inferStrain(getBetas(sm, mask=FALSE))$probs, decreasing=TRUE)
predictAgeHorvath353(getBetas(sh, mask=FALSE))
predictMouseAgeInMonths(getBetas(sm, mask=FALSE))
predictMouseAgeInMonths(getBetas(dyeBiasNL(noob(sm)), mask=FALSE))
sh.normal = sesameDataGet(“EPIC.5.SigDFs.normal”)
length(sh.normal)
segs = cnSegmentation(sh, sh.normal)
visualizeSegments(segs)
betas = sesameDataGet(“HM450.1.TCGA.PAAD”)$betas
estimateLeukocyte(betas)
compareMouseTissueReference(getBetas(sm))

## For more information, visit https://www.bioconductor.org/packages/release/bioc/html/sesame.html
