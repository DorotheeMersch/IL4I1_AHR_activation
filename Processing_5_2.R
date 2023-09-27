library(oligo)
library(AffyGEx)
library(tidyverse)
library(limma)
source('/Users/dorimersch/Desktop/Dorothee/functions_processing.R') 
source("/Users/dorimersch/AHR_IL4I1_workflows/AHR_IL4I1_manuscript/functions_DG.R")

AHR_genes <- read.delim("/Users/dorimersch/AHR_IL4I1_workflows/AHR_IL4I1_manuscript/Resources/overlapping_AHR_signature_genes.txt", sep = "\t")

# read cel files
setwd('/Users/dorimersch/Desktop/Dorothee/5_affy_chips01')
celfiles <- list.files(pattern = "\\.CEL$")
raws_cel <- read.celfiles(celfiles)

# read phenoData
pdat <- read_excel("~/Desktop/Dorothee/RNA_for_HD_markiert_Saskia_sorted_by_color.xlsx", sheet = 6)
                                           
# adjust column names, specify date and experimentator (lysis, isolation etc.)
## add ondition column to pdat
pdat$condition <- gsub("^.*\\_","",pdat$LINA) %>% gsub(" ", "_",.)%>%paste0("cond_",.)
pdat$donor <- gsub("_CD4.*","",pdat$LINA) %>% ifelse(nchar(.)==0,"unknown",.)%>%gsub("_Spender ","",.)

## write phenoData only with information needed (date, sampleID, condition)
pdat_filtered <- data.frame('date_1'=pdat$Datum_Ernte, 'sample_id'=pdat$`Probe-ID`, 'condition'=pdat$condition, 
                            'date_2'=pdat$Datum_2, 'experimentator'=pdat$Experimentator_2, 'donor'=pdat$donor)
rownames(pdat_filtered) <- pdat$`Probe-ID`

## get phenoData in same order as raws_cel 
sampleNames(raws_cel) <- sub("_.*", "", sampleNames(raws_cel))
pdat_filtered <- pdat_filtered[match(sampleNames(raws_cel), rownames(pdat_filtered)),]

## create phenoData AnnotatedDataFrame
metadata <- data.frame(labelDescription=c('date_1', 'sample_id', 'condition', 'date_2', 'experimentator', 'donor'), channel=factor('_ALL_'))
pd <- new('AnnotatedDataFrame', data=pdat_filtered, varMetadata=metadata)
## define phenoData in raws_cel 
phenoData(raws_cel) <- pd


# use gen_rma_FUN from AffyGEx package to:
## - convert affyBatch into ExpressionSet
## - generate featureData
## - assign gene symbols
## - reduce probesets to genesets
## - create gene set values
dir.create("./RDS")
dir.create("./TopTables")
dir.create("./GSA")
rma_2sd <- gen_rma_FUN(raws_cel, sd_val = 2,snames = sampleNames(raws_cel),ct_off = 0.75)

## PCA before correction
PCA_raw <- prcomp(t(exprs(rma(raws_cel))))

df_PCA_raw <- data.frame(PCA_raw$x[,1:10],pd@data)
df_PCA_raw$date_1 <- as.factor(df_PCA_raw$date_1)
df_PCA_raw$date_2 <- as.factor(df_PCA_raw$date_2)
list_PCs <- list(c('PC1','PC2'),c('PC1','PC3'), c('PC1','PC4'), c('PC2','PC3'), c('PC2','PC4'), c('PC3','PC4'))

pdf(paste("./Figures/EXP_5_FIGURE_PCA_no_correction.pdf",sep = ""),width = 12, height = 8)
for (variable in colnames(pdat_filtered)[-2]) {
  for (PCs in list_PCs) {
    print(ggscatter(df_PCA_raw, PCs[1],PCs[2],color=variable)+theme_pubr())
  }
}
dev.off()

## PCA after correction
PCA_rma <- prcomp(t(exprs(rma_2sd)))

df_PCA_rma <- data.frame(PCA_rma$x[,1:10],pd@data)
df_PCA_rma$date_1 <- as.factor(df_PCA_rma$date_1)
df_PCA_rma$date_2 <- as.factor(df_PCA_rma$date_2)

pdf(paste("./Figures/EXP_5_FIGURE_PCA_after_correction_rma.pdf",sep = ""),width = 12, height = 8)
for (variable in colnames(pdat_filtered)[-2]) {
  for (PCs in list_PCs) {
    print(ggscatter(df_PCA_rma, PCs[1],PCs[2],color=variable)+theme_pubr())
  }
}
dev.off()


# use updated anova function to compute pov (from functions_processing.R)
PoV_PCs_list = list()
for (variable in colnames(pdat_filtered)[-2]) {
  PoV_PCs_list[[variable]] = pca_anova(rma_2sd, variable)
}
PoV_PCs_all = do.call(rbind, PoV_PCs_list)


# MAYBE CHECK INDEPENDENCE OF VARIABLES VIA CONTINGENCY TABLE AND CHI SQAURE TEST
# chisq_pval is a function itegrated in functions_processing.R
## for all combinations of variables
chisq_list1 = list()
chisq_list2 = list()
for (variable1 in colnames(pdat_filtered[-2])) {
  for(variable2 in colnames(pdat_filtered[-2])) {
    chisq_res = chisq_pval(rma_2sd, variable1, variable2)
    chisq_list1[[variable2]] = chisq_res
  } 
  #output[variable1,variable2]=chisq_res
  chisq_list2[[variable1]]=chisq_list1
}

## append PoV_PCs_all df with chisq results
chisq_df <- as.data.frame(do.call(rbind, chisq_list2))
POV_PCs_chisq <- cbind(PoV_PCs_all, chisq_df)


#DGE
# Design matrix
## choose one of date_1, date_2 and experimentator -> they are all dependent from each other (see chisq results)
design_eset <- model.matrix(~0+condition+date_1, data=pd@data)
design_eset <- as.matrix(design_eset)

## block donor
corfit <- duplicateCorrelation(rma_2sd,design_eset,block=phenoData(rma_2sd)$donor)

arr_wts <- arrayWeights(rma_2sd,design = design_eset)

colnames(design_eset) <- gsub("condition","",colnames(design_eset))
contrast_matrix <- makeContrasts(I3P=cond_3uM_I3P-cond_005_DMSO, HPP=cond_40uM_HPP-cond_005_DMSO,
                                 levels = design_eset)
## DGE 2sd
fit_2sd <- lmFit(rma_2sd, design = design_eset, block = phenoData(rma_2sd)$donor, correlation = corfit$cor, weights = arr_wts)
fit_2sd_cont <- contrasts.fit(fit_2sd, contrasts = contrast_matrix)
fit_2sd_eB <- eBayes(fit_2sd_cont, trend = T)
saveRDS(fit_2sd_eB, "./RDS/fit_2sd_eB_EXP5.rds")

## all tts
tts_2sd_u <- map(1:2, function(i){
  tt <- topTable(fit_2sd_eB, coef = i,  sort.by = "M",  number = Inf, adjust.method = "BH")
  tt <- tt[,c(29,25,40:45)]
})
names(tts_2sd_u) <- colnames(contrast_matrix)

## annotation update
library(dplyr)
tts_2sd_u[[1]] <- annot_hgcn_FUN2(tts_2sd_u[[1]], hgnc = hgnc) #annot_hgcn_FUN2 from functions_DG.R
tts_2sd_u[[2]] <- annot_hgcn_FUN2(tts_2sd_u[[2]], hgnc = hgnc)

saveRDS(tts_2sd_u, "./RDS/tts_2sd_u_annotated.rds")

## add AHR signature genes
tts_2sd_u_AHR <- map(tts_2sd_u, function(x,y){
  x$AHR_target <- match(x$hGene, y)
  x$AHR_target[!is.na(x$AHR_target)] <- "yes"
  x
}, y=AHR_genes$Gene)

walk2(tts_2sd_u_AHR, paste("./TopTables/tt_2sd_",names(tts_2sd_u),"_AHR.txt",sep = ""), write.table, quote = F, sep = "\t")

## Barcodeplots AHR signature
index_1 <- tts_2sd_u_AHR$I3P$hGene%in%AHR_genes$Gene
index_2 <- tts_2sd_u_AHR$HPP$hGene%in%AHR_genes$Gene

pdf("./Figures/barcodeplot_EXP_5.pdf")
barcodeplot(tts_2sd_u_AHR$I3P$t,index_1, main="3uM I3P in CD4 PBMC")
barcodeplot(tts_2sd_u_AHR$HPP$t,index_2, main="40 uM HPP in CD4 PBMC")
dev.off()


#### GSA ####
design_camera <-  model.matrix(~0+condition+date_1+donor, data=pd@data)
## Roast for the AHR signature
roast_FUN <- function(tt_ids, cont_mat,gs,rma_obj, des_mat, arr_wts=NULL, ...){
  AHR_indexed_hi_lo <- ids2indices(gene.sets = gs, identifiers = tt_ids[,1])
  prbsets_gsa_2sd <- rownames(tt_ids)
  gsa_2sd_idx <- match(prbsets_gsa_2sd, rma_obj@featureData@data$probesetid)
  gsa_2sd_eset <- rma_obj[gsa_2sd_idx]
  AHR_hi_lo_roast <- roast(oligo::exprs(gsa_2sd_eset),weights=arr_wts,index = AHR_indexed_hi_lo, design = des_mat,
                           contrast = cont_mat, set.statistic = "floormean")
  AHR_hi_lo_roast
}

all_roast <- map2(tts_2sd_u_AHR, as.data.frame(contrast_matrix), roast_FUN, gs=AHR_genes$Gene, rma_obj=rma_2sd, des_mat=design_eset, arr_wts=arr_wts, block = donor, correlation = corfit$cor) %>% do.call(rbind,.)
all_roast <- data.frame(Condition=rownames(all_roast), all_roast, stringsAsFactors = F)

write.table(all_roast, "./GSA/AHR_enrichment_I3P_HPP_all.txt", row.names = F, sep = "\t")

gsa_res <- gsa_FUN2(top_t = tts_2sd_u_AHR$I3P, rma_obj = rma_2sd, exprmnt = "Exp_5", cntrst = "many",
                   cnt_mat = contrast_matrix, res_path = "./GSA/",wts = arr_wts,
                   msig.data.lists = msig.data.lists, msigs = 14, pltfrm = "Affy", d_m = design_eset)

for(i in 1:length(gsa_res)){
  walk(1:length(gsa_res[[i]]), safely(gsa_bar_plot_FUN), gsa_ls=gsa_res[[i]], res_path="/Figures/",
       coi=names(gsa_res)[i], wd=12, ht=8, file_type="pdf")
}

gsa_paths <- list.files("./GSA/",pattern = "Exp_5",full.names = T)
gsa_paths <- gsa_paths[grep("\\.txt",gsa_paths)] %>% .[-1]
All_tts_GSAs <- map(gsa_paths,read.delim,stringsAsFactors=F)
gsa_paths_names <- gsub("\\.//|\\.txt","",gsa_paths)
names(All_tts_GSAs) <- gsa_paths_names

condition_names<-gsub("gsa_|_Exp_5.*.","",gsa_paths_names)

All_tts_GSAs_05 <- map2(All_tts_GSAs,condition_names, function(x,y){
  df <- data.frame(condition=y, pathways=rownames(x),x,stringsAsFactors = F)
  df <- df[df$PValue<=0.05,]
  df$Direction <- factor(df$Direction, levels = c("Up","Down"))
  df$logPV <- -log10(df$PValue)
  if(length(df$pathways)>30){
    df[1:15,]
  } else {
    df
  }
})

names(All_tts_GSAs_05) <- gsub("\\./GSA//","",names(All_tts_GSAs_05))


## Volcano plots
tt_I3P<- tt_AHR_FUN(tts_2sd_u_AHR$I3P, AHR_genes$Gene, 0.58, 2)
v_p_I3P <- volc_plot(tt_I3P)

pdf("./Figures/EXP5_I3P_volcano_plot.pdf", width = 8, height = 8, pointsize = 28)
v_p_I3P
dev.off()

tt_HPP<- tt_AHR_FUN(tts_2sd_u_AHR$HPP, AHR_genes$Gene, 0.58, 2)
v_p_HPP <- volc_plot(tt_HPP)

pdf("./Figures/EXP5_HPP_volcano_plot.pdf", width = 8, height = 8, pointsize = 28)
v_p_HPP
dev.off()