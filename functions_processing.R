# functions for pre-processing
library(sjstats)
library(sva)

# create pca plot
pca_colored_by_variable <- function(exprs_matrix, variable, PC_x, PC_y){
  # perform pca
  PCA <- prcomp(exprs_matrix, scale = FALSE)
  
  ## plot PCA:
  ### - colored by condition
  ### - colored by date
  variable_for_color = factor(variable)
  
  plot_PCA_colored_condition <- ggplot(as.data.frame(PCA$rotation), aes(x= PCA$rotation[,PC_x], y = PCA$rotation[,PC_y]))+
    geom_point(aes(colour=variable_for_color))

  return(plot_PCA_colored_condition)
}

# compute PoV
PoV_explained_anova <- function(rma_object, variable_name){
  pca_new <- prcomp(t(exprs(rma_object)[,rownames(rma_object@phenoData@data)]))
  
  pca_anova = apply(t(pca_new$x), 1, function(x, targets){
    PCs = x
    factor_ids = factor(targets[[variable_name]])
    anova = aov(PCs ~ factor_ids)
    aov_stats = anova_stats(anova)
  },targets = phenoData(rma_object))
  
  pca_anova <- pca_anova %>% bind_rows(.id = "PC") %>% as_tibble() %>% 
    group_by(PC) %>% gather(stats,value,-(PC:term))  %>% spread(term,value) %>%
    ungroup()
  
  pca_res = summary(pca_new)
  PCs_date = pca_anova %>% filter(stats == "p.value" & factor_ids < 0.05)
  
  PoV <- as.data.frame(pca_res$sdev^2 / sum(pca_res$sdev^2),
                       row.names = colnames(pca_res$rotation))
  sum_PoV_PCs <- sum(subset(PoV, rownames(PoV) %in% PCs_date$PC))
  sum_PoV_PCs
}


PoV_explained_ttest <- function(rma_object, variable, value_1, value_2){
  PCA <- prcomp(exprs(rma_object), scale = FALSE)
  
  rma_object@phenoData@data[[variable]] <- as.character(rma_object@phenoData@data[[variable]])
  
  samples_1 <- as.data.frame(dplyr::filter(phenoData(rma_object)@data, rma_object@phenoData@data[[variable]] == value_1))
  PCs_1 <- subset(PCA$rotation, rownames(PCA$rotation) %in% rownames(samples_1))
  
  samples_2 <- as.data.frame(dplyr::filter(phenoData(rma_object)@data, rma_object@phenoData@data[[variable]] == value_2))
  PCs_2 <- subset(PCA$rotation, rownames(PCA$rotation) %in% rownames(samples_2))
  
  pvalues <- sapply(colnames(PCA$rotation), function(x){
    t.test(PCs_1[,x], PCs_2[,x])$p.value
  })
  
  PCs_variable <- rownames(as.data.frame(pvalues) %>%
                             rownames_to_column('PCs') %>%
                             filter(pvalues < 0.05) %>%
                             column_to_rownames('PCs'))
  
  PoV <- as.data.frame(PCA$sdev^2 / sum(PCA$sdev^2), row.names = colnames(PCA$rotation))
  sum_PoV <- sum(subset(PoV, rownames(PoV) %in% PCs_variable))
  return(sum_PoV)
}


# new function adapted from PoV_explained_anova 
pca_anova <- function(rma_object, variable){
  pca_new <- prcomp(t(exprs(rma_object)[,rownames(rma_object@phenoData@data)]))
  pca_res <- summary(pca_new)
  factor_ids = factor(phenoData(rma_object)[[variable]])
  anova = aov(pca_new$x ~ factor_ids)
  anova_summary = summary(anova)
  df_anova <- do.call("rbind", anova_summary) %>% na.omit(.) 
  rownames(df_anova) <- colnames(pca_new$x)
  df_anova$PoV <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  PoV_PCs <- data.frame('PoV_sum' = sum(df_anova[df_anova$`Pr(>F)` < 0.05,]$PoV), 
                        'PCs_sum' = toString(rownames(df_anova[df_anova$`Pr(>F)` < 0.05,])))
  rownames(PoV_PCs) = variable
  return(PoV_PCs)
}

# function to test dependence between categorical variables using chisquare test
chisq_pval <- function(rma_object, variable1, variable2){
  cont_table = table(phenoData(rma_object)@data[,c(variable1, variable2)])
  pval = chisq.test(cont_table, correct=F)$p.value
  return(pval)
}