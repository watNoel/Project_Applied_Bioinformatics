import numpy as np
import pandas as pd

go_data = pd.read_csv('Project_Applied_Bioinformatics/GO_terms_association.txt',delimiter= '\t', names=range(11))

col_names = ['GO function', 'gene']
functions = []
genes_grouped_by_function = []
for row in range(0, len(go_data[0])):
  functions.append(go_data[1][row])
  genes=[]
  
  for gene in go_data[row][10]:
    gene_list = gene.split(',')

  genes_grouped_by_function.append(gene_list)
  
    
    
    
  

  
  
  
  
