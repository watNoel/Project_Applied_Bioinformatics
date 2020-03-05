
library(ggrepel)
library(EnhancedVolcano)
library(farver)
library(forcats)
library(tidyverse)

# The required files are 3, a file with the GCN4 genes in a newline separated list,
#and the output file cds_exp.diff from cuffdiff, 
# and a file with associated GO-terms. Note that the script is highly specific to the format of these files
# due to the manipulation of the tables to combine the GO-term association with the differential-expression analysis
# for example, note that the Go-term association file is assumed to be separated by ", ", a comma+space 



###-----  TO GET THE TWO SEPARATE PLOTS, Change the filename in read.csv*       ###
###-----  AND change the name of the TITLE in the arguments to EnhancedVolcano. ###

#diff_expr <- read.csv(file="cds_exp_diff_wt", sep= '\t') for wildtype vs wildtype
diff_expr <- read.csv(file="exp_diff_genes_wt_vs_deletion", sep= '\t')

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


# this is some shady loop to introuce the "everything else label"
levels(diff_expr_and_goterms$GOterm)=c(levels(diff_expr_and_goterms$GOterm),"Everything else")

for (i in 1:length(diff_expr_and_goterms$gene)){
  if (is.na(diff_expr_and_goterms$GOterm[i]) ) {
    diff_expr_and_goterms$GOterm[i] = as.factor("Everything else") 
  }
}



#### ------------ This segment gives the custom categories of genes
keyvals= ifelse(diff_expr_and_goterms$GOterm =="fungal-type cell wall organization or biogenesis",
                'darkgreen',ifelse(diff_expr_and_goterms$GOterm == "cellular amino acid metabolic process",
                       'firebrick3', ifelse(diff_expr_and_goterms$GOterm== "glycolytic process",
                                       'darkmagenta','gray87')
                              )
               )  
                  
names(keyvals)[keyvals == 'gray87'] <- 'everything else'
names(keyvals)[keyvals == 'darkgreen'] <- 'cell wall organization or biogenesis'
names(keyvals)[keyvals == 'firebrick3'] <- 'cellular amino acid metabolic process'
names(keyvals)[keyvals == 'darkmagenta'] <- 'glycolytic process'



###3 --- this segment gives the cutom shape for specific GCN4-target genes
# requires input file which ONLY consists of a list of the relevant gene-names. 
# this list was aquired by a´taking all the common entries in the lists of genes from cuffdiff
# and https://www.yeastgenome.org/locus/GCN4/regulation under gene targets.  

GCN4_genes <- read_lines("./GCN4/GCN4_targets_in_diff_exp_output")

keyvals.shape <- ifelse( diff_expr_and_goterms$gene %in% GCN4_genes, 17,19) # is triangle

names(keyvals.shape)[keyvals.shape==17] = 'GCN4-associated genes'
names(keyvals.shape)[keyvals.shape==19] = ''

            


####--------- Finally, I add 1 to each value and calculate log-foldchanges to avoid infinite changes   

diff_expr_and_goterms$value_1=diff_expr_and_goterms$value_1+1
diff_expr_and_goterms$value_2=diff_expr_and_goterms$value_2+1
diff_expr_and_goterms$adj_log2fold= log2(diff_expr_and_goterms$value_2/diff_expr_and_goterms$value_1)




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
                shapeCustom = keyvals.shape,
                #title = 'WT  (1.3% Isobutanol)  compared to WT (0% Isobutanol)', # for WtvsWT
                title = 'DGLN3  (1.3% Isobutanol)  compared to WT  (1.3% Isobutanol)',
                legendPosition = 'right',
                legendIconSize = 5.0,
                legendLabSize = 8,
                axisLabSize = 14
)


###







# From here it does not work, old code ----------------------------
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

