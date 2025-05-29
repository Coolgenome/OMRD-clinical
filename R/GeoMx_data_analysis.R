### GeoMx data analysis ###

### load essential packages ###
library(tidyverse)
library(dplyr)
library(ggforce)
library(ggplot2)
library(ggalluvial)
library(ggrepel)
library(reshape2)
library(cowplot)
library(fgsea)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(NanoStringNCTools)
library(SpatialOmicsOverlay)
library(GeomxTools)
library(GeoMxWorkflows)

### Data QC followed standard pipelines by Nanostring, with the package GeomxTools ###

#DEG analysis
# run LMM:
# formula follows conventions defined by the lme4 package
library("optimx")
ind <- pData(MRD_Data)$SLL_outcome == "Positive" & pData(MRD_Data)$SegmentLabel == "Tumor" & pData(MRD_Data)$Timepoint %in% c("SLL","Interval_debulking")

mixedOutmc <- mixedModelDE(MRD_Data[, ind],
                           elt = "log_q",
                           modelFormula = ~ Timepoint + (1 | SlideName),
                           groupVar = "Timepoint")

# format results as data.frame
r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
tests <- rownames(r_test)
r_test <- as.data.frame(r_test)
r_test$Contrast <- tests

# use lapply in case you have multiple levels of your test factor to correctly associate gene name with it's row in the results table
r_test$Gene <- unlist(lapply(colnames(mixedOutmc), rep, nrow(mixedOutmc["lsmeans", ][[1]])))
r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
r_test <- r_test[, c("Gene", "Contrast", "Estimate", "Pr(>|t|)", "FDR")]
results <- r_test

# Categorize Results based on P-value & logFC for plotting
results$Color <- "NS"
results$Color[results4$FDR < 0.1 & abs(results4$Estimate) >= 0.5] <- "FDR<0.1 & |logFC>=0.5|"
results$invert_P <- (-log10(results$FDR)) 

# Volcano plot, with part of the DEGs labeled
ggplot(results, aes(x = Estimate, y = invert_P, color = Color,label = Gene)) +
  geom_point(size=0.3) +
  labs(x = "Enriched in Interval_debulking <- log2(FC) -> Enriched in SLL",
       y = "-log10(p.adj)") +
  scale_color_manual(values = c(
    `FDR < 0.1` = "#e31a1c",
    `NS` = "#d9d9d9"),
    guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, (FDR<0.05 & abs(Estimate)>=1) | Gene %in% c("DUSP1","GPI","KDM3A")),
                  size = 3, point.padding = 0.2, color = "black",
                  min.segment.length = .5, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") 

#GSEA
## prepare ranked gene list
SLL_Pos_SLL_vs_Interval_genes <- results %>% arrange(desc(Estimate)) %>% dplyr::select(Gene,Estimate)
SLL_Pos_SLL_vs_Interval_genes <- tibble::deframe(SLL_Pos_SLL_vs_Interval_genes)

## prepare geneset
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
hs_msigdb_df_filter <- filter(hs_msigdb_df, hs_msigdb_df$gs_cat=="H") %>% dplyr::select("gs_name","human_gene_symbol")
msigdb_hallmark <- split(hs_msigdb_df_filter$human_gene_symbol, hs_msigdb_df_filter2$gs_name)

for (i in 1:length(msigdb_hallmark)) {
  geneset <- msigdb_hallmark[[i]]
  geneset <- unique(geneset)
  msigdb_hallmark[[i]] <- geneset
}


## perform GSEA
set.seed(42)
SLL_Pos_SLL_vs_Interval_fgseaRes <- fgsea(pathways=msigdb_hallmark, stats=SLL_Pos_SLL_vs_Interval_genes, nPermSimple = 10000)

SLL_Pos_SLL_vs_Interval_sig_pathway <- SLL_Pos_SLL_vs_Interval_fgseaRes[padj < 0.05][head(order(NES, decreasing = T), n=50), pathway]

plotGseaTable(msigdb_hallmark[SLL_Pos_SLL_vs_Interval_sig_pathway], SLL_Pos_SLL_vs_Interval_genes, SLL_Pos_SLL_vs_Interval_fgseaRes,gseaParam=0.5)

plotEnrichment(msigdb_hallmark[["HALLMARK_HYPOXIA"]], SLL_Pos_SLL_vs_Interval_genes) + 
  labs(title="Hypoxia")




### analysis of trend of gene expression along timepoints ###
exp_mat<-as.data.frame(t(MRD_Data@assayData$log_q))
exp_mat<-rownames_to_column(exp_mat,var="AOI_barcode")
annotation<-rownames_to_column(annotation,var="AOI_barcode")
exp_mat<-left_join(exp_mat, annotation, by="AOI_barcode")
exp_mat_SLL_pos_tumor<-filter(exp_mat_SLL_pos, exp_mat_SLL_pos$SLL_outcome=="Positive" & exp_mat_SLL_pos&SegmentLabel=="Tumor")


### plot gene expression across timepoints ###
ggplot(exp_mat_SLL_pos_tumor, aes(x=Timepoint, y= ABCA1))+
  geom_violin(aes(fill=Timepoint), scale = "width",adjust=1)+
  geom_boxplot(width=0.3,outlier.size = 0)+
  geom_point(color="lightgrey",size=1.5,position="jitter")+
  scale_fill_manual(values = c("#9ecae1","#4292c6","#08306b"))+
  scale_color_manual(values = c("#9ecae1","#4292c6","#08306b"))+
  theme_classic()+
  theme(plot.title = element_text(size=12), 
        axis.title = element_text(size=12), 
        axis.text.x = element_text(angle = 45, hjust=1), 
        legend.position = "none")
  ylab("Log normalized gene expression")


## Baysian model for estimating the monotonic trend of gene expression  ##  
library(brms)
  
# scale expression data, and add metadata #
exp_mat_SLL_pos_tumor_exp <- exp_mat_SLL_pos_tumor[,2:14411]
rownames(exp_mat_SLL_pos_tumor_exp) <- exp_mat_SLL_pos_tumor[,1]
exp_mat_SLL_pos_tumor_exp_scale <- scale(exp_mat_SLL_pos_tumor_exp)
exp_mat_SLL_pos_tumor_exp_scale <- cbind(exp_mat_SLL_pos_tumor_exp_scale, exp_mat_SLL_pos_tumor[,14412:14503])
exp_mat_SLL_pos_tumor_exp_scale <- rownames_to_column(exp_mat_SLL_pos_tumor_exp_scale, var="AOI_barcode")
exp_mat_SLL_pos_tumor_exp_scale$Timepoint <- factor(exp_mat_SLL_pos_tumor_exp_scale$Timepoint, 
                                                    levels = c("Pre_NACT", "Interval_debulking", "SLL"),
                                                    ordered = TRUE)

brm_model_res<-data.frame(Gene=rownames(target_MRD_Data_filter2@assayData$log_q),estimate=NA, stderr=NA, ci_lower=NA, ci_upper=NA)

for (i in 1:nrow(brm_model_res)) {
  gene_to_test = brm_model_res$Gene[i]
  df <- exp_mat_SLL_pos_tumor_exp_scale[,c(gene_to_test,"AOI_barcode","Timepoint","SlideName")]
  colnames(df)[1] <- "expression"
  model<-brm(expression ~ mo(Timepoint) + (1 | SlideName), data = df)
  summary_df <- as.data.frame(fixef(model))
  mo_row <- summary_df[rownames(summary_df) == "moTimepoint", ]
  brm_model_res[i,"estimate"] = mo_row$Estimate
  brm_model_res[i,"stderr"] = mo_row$Est.Error
  brm_model_res[i, "ci_lower"] = mo_row$Q2.5
  brm_model_res[i, "ci_upper"] = mo_row$Q97.5
}



