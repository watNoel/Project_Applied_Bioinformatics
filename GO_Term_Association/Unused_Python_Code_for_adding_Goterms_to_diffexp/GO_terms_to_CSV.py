import numpy as np
import pandas as pd

def go_to_df(file_path_input, file_path_output = None):
    """
    Takes the Gene Ontolpgy (GO) annotation file and creates a dataframe with all genes aswell as their GO annotation.
    Optional: Write a CSV file
    """
    go_data = pd.read_csv(file_path_input,delimiter='\t')

    functions = []
    genes_grouped_by_function = []

    for row in range(go_data.shape[1]):
        functions.append(go_data["TERM"][row])
        gene_list = go_data['ANNOTATED_GENES'][row].split(',')
        genes_grouped_by_function.append(gene_list)

    ordered_terms = []
    ordered_genes = []
    for term in range(len(functions)):
        for gene in range(len(genes_grouped_by_function[term])):
            ordered_terms.append(functions[term])
            gene_no_whitespace = genes_grouped_by_function[term][gene]
            ordered_genes.append(gene_no_whitespace.strip())

    go_term_gene = pd.DataFrame({"GO term":ordered_terms, "gene":ordered_genes})       

    if file_path_output != None:

        go_term_gene.to_csv(file_path_output)
    
    return go_term_gene

go_term_gene = go_to_df('GO_terms_association.txt', 'go_terms_and_gene_names.csv')

# Read the cuffdiff output (wt and gln3d) file and merge with gene annotation file
gene_expression_wt = pd.read_csv('cuffdiff_output/cdout/wt/gene_exp.diff', delimiter= '\t')
merged_df_wt = gene_expression_wt.merge(go_term_gene, on='gene', how='left')
merged_df_wt.to_csv('wt_diff_gene_expression_with_GO.csv')

gene_expression_wt_dgln3 = pd.read_csv('cuffdiff_output/cdout/wt_dgln3/gene_exp.diff', delimiter= '\t')
merged_df_wt_dgln3 = gene_expression_wt_dgln3.merge(go_term_gene, on='gene', how='left')
merged_df_wt_dgln3.to_csv('dgln3_wt_diff_gene_expression_with_GO.csv')
