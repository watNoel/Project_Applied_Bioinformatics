#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')
#
#BiocManager::install('EnhancedVolcano')

#install.packages("farver")
library(ggrepel)
library(EnhancedVolcano)
library(farver)
## begin by importing the file into R as a data.frame

diff_expr <- read.csv(file="cds_exp_diff_wt", sep= '\t')

# we only care about a few columns of this file! 

frame_for_volcano <- diff_expr[,c("gene","value_1","value_2","log2.fold_change.","p_value","q_value")]

frame_for_volcano$value_1=frame_for_volcano$value_1+1
frame_for_volcano$value_2=frame_for_volcano$value_2+1
frame_for_volcano$adj_log2fold= log2(frame_for_volcano$value_2/frame_for_volcano$value_1)
# now we have all things nicely formattted
rownames(frame_for_volcano)= diff_expr$gene
#[,c("adj_log2fold","p_value"]

EnhancedVolcano(frame_for_volcano,
                lab = rownames(frame_for_volcano),
                x = 'adj_log2fold',
                y = 'p_value',
                xlim = c(-10, 10),
                ylim = c(0,5.5),
                pCutoff = 0.05,
                FCcutoff = 1)

frame_for_volcano$GO_term<-0

GOterm_enrichment <- read.csv("GO_terms", sep="\t")


Go_terms <- read.table(file="only_genes.txt")

rownames(Go_terms)= GOterm_enrichment$TERM


