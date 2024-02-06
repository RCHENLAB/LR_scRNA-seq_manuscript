import pandas as pd

# Read in the cell class by isoform count matrix and SQANTIs isoform classification
list1 = pd.read_csv('./LR.flair.quantify', sep = ",", encoding = "ISO-8859-1", header = (0))
list2 = pd.read_csv('./SQANTI3_classification.txt', sep = "\t", encoding = "ISO-8859-1", header = (0))
df = list1.merge(list2, how='inner', on=["id"])
df['sum'] = df['AC'] + df['BC'] + df['Cone'] + df['MG'] + df['RGC'] + df['Rod']

# Filter out isoforms with no more than 5 UMIs supported and the ones classified as artifacts by SQANTI3
df = df.loc[(df['sum'] > 5)] & df.loc[(df['filter_result'] == 'Isoform')]

# Filter out isoforms from genes with a UMI count of 10 or less
df1 = df.groupby(['gene']).size().reset_index(name='count')
df2 = df.groupby('gene')['sum'].sum()
isoform_counts = pd.merge(df1, df2, on=['gene'])
isoform_counts = isoform_counts.loc[(isoform_counts['sum'] > 10)]
df3 = pd.merge(df, isoform_counts[['gene', 'sum', 'count']], on=['gene'])

# Calculate the poportion of isoform abundance in the gene and filter out isoforms with a UMI count less than 2% of the total UMI count for the respective gene
df3['iso_exp_percentage_gene'] = df3['sum_x'] / df3['sum_y']
df3 = df3.loc[(df3['iso_exp_percentage_gene'] > 0.02)]

# Re-calculate the poportion of isoform abundance in the gene after filtering
df4 = df3.groupby(['gene']).size().reset_index(name='count')
df5 = df3.groupby('gene')['sum'].sum()
isoform_counts2 = pd.merge(df4, df5, on=['gene'])
final = pd.merge(df3, isoform_counts2[['gene', 'sum', 'count']], on=['gene'])
final['iso_exp_percentage_gene_after'] = final['sum_x'] / final['sum_y']

# Write the output to a file
final.to_csv("./LR.flair.quantify_groupGene.csv", index = False, sep = ",")
