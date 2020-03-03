#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')
#
#BiocManager::install('EnhancedVolcano')

#install.packages("farver")
library(ggrepel)
library(EnhancedVolcano)
library(farver)
library(forcats)
## begin by importing the file into R as a data.frame



# experimental stuff, but it works! ------------------------------------------------------




library(tidyverse)

diff_expr <- read.csv(file="cds_exp_diff_wt", sep= '\t')
GOterm_enrichment <- read.csv("association_for_volcano_plot.txt", sep="\t", header= 0) %>%
  as.tibble() %>% 
  mutate(gene= strsplit(as.character(V11),", ")) %>% 
  unnest(gene)

diff_expr_and_goterms <- full_join(
GOterm_enrichment %>% 
  distinct(V2,gene),
diff_expr %>% mutate( gene = as.character(gene))
, by = "gene" 
  ) %>% rename(GOterm=V2)

levels(diff_expr_and_goterms$GOterm)=c(levels(diff_expr_and_goterms$GOterm),"Everything else")

for (i in 1:length(diff_expr_and_goterms$gene)){
  if (is.na(diff_expr_and_goterms$GOterm[i]) ) {
    diff_expr_and_goterms$GOterm[i] = as.factor("Everything else") 
  }
}


keyvals= ifelse(diff_expr_and_goterms$GOterm =="fungal-type cell wall organization or biogenesis",
                'green',ifelse(diff_expr_and_goterms$GOterm == "cellular amino acid metabolic process",
                       'red', ifelse(diff_expr_and_goterms$GOterm== "glycolytic process",
                                       'blue','grey')
                              )
               )  
                  
           

names(keyvals)[keyvals == 'grey'] <- 'everything else'
names(keyvals)[keyvals == 'green'] <- 'fungal-type cell wall organization or biogenesis'
names(keyvals)[keyvals == 'red'] <- 'cellular amino acid metabolic process'
names(keyvals)[keyvals == 'blue'] <- 'glycolytic process'
               
diff_expr_and_goterms$value_1=diff_expr_and_goterms$value_1+1
diff_expr_and_goterms$value_2=diff_expr_and_goterms$value_2+1
diff_expr_and_goterms$adj_log2fold= log2(diff_expr_and_goterms$value_2/diff_expr_and_goterms$value_1)


GCN4_genes <- read_lines("./GCN4/GCN4_targets_in_diff_exp_output")

keyvals.shape <- ifelse(diff_expr_and_goterms$gene %in% GCN4_genes, 17, 1)

names(keyvals.shape)[keyvals.shape==17] = 'GCN4-associated genes'
names(keyvals.shape)[keyvals.shape==1] = ''

EnhancedVolcano(diff_expr_and_goterms,
                lab = diff_expr_and_goterms$gene,
                x = 'adj_log2fold',
                y = 'p_value',
                selectLab = "", # want no labels!
                xlim = c(-10, 10),
                ylim = c(0,5.5),
                # xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 1,
                colCustom = keyvals,
                shapeCustom = keyvals.shape
)





# From here it does not seem to work, old code ----------------------------





#GOterm_enrichment_tsv <- read_tsv("association_for_volcano_plot.txt", colnames= FALS)



#Go_terms <- read.table(file="only_genes.txt",header = 0)


diff_expr <- read.csv(file="wt_diff_gene_expression_with_GO.csv", sep= ',')
diff_expr <- read.csv(file="cds_exp_diff_wt", sep= '\t')

# we only care about a few columns of this file! 

frame_for_volcano <- diff_expr[,c("gene","value_1","value_2","log2.fold_change.","p_value","q_value", "GO.term")]
levels(frame_for_volcano$GO.term)=c(levels(frame_for_volcano$GO.term),"Everything else")

# this creates the last factor
for (i in 1:length(frame_for_volcano$GO.term)){
  if (frame_for_volcano$GO.term[i] == "") {
    print("hi")
    frame_for_volcano$GO.term[i] = as.factor("Everything else") 
  }
  
}

#frame_for_volcano$GO.term=forcats::fct_explicit_na(frame_for_volcano$GO.term)

frame_for_volcano$value_1=frame_for_volcano$value_1+1
frame_for_volcano$value_2=frame_for_volcano$value_2+1
frame_for_volcano$adj_log2fold= log2(frame_for_volcano$value_2/frame_for_volcano$value_1)
# now we have all things nicely formattted
#rownames(frame_for_volcano)= diff_expr$gene
#[,c("adj_log2fold","p_value"]

keyvals <- ifelse(
  frame_for_volcano$adj_log2fold < -2.5, 'royalblue',
  ifelse(frame_for_volcano$adj_log2fold > 2.5, 'gold',
         'black'))
names(keyvals)[keyvals == 'gold'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'



EnhancedVolcano(frame_for_volcano,
                lab = frame_for_volcano$gene,
                x = 'adj_log2fold',
                y = 'p_value',
                selectLab = "", # want no labels!
                xlim = c(-10, 10),
                ylim = c(0,5.5),
                # xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 1,
                colCustom = keyvals,
)


keyvals2=as.character(frame_for_volcano$GO.term)
names(keyvals2)=as.character(frame_for_volcano$GO.term)

