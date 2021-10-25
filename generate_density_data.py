# import libraries
import sys
import io
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Argvariables
#print('arguments', sys.argv)
vcf_file = str(sys.argv[1])
window_size = int(sys.argv[2])
increment_value = int(sys.argv[3])


# Read vcf file without headers ==============================================
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    

df = read_vcf(vcf_file)

# cleaning data ==============================================
## format CHROM column
df['CHROM'] = df['CHROM'].str.replace('chr0','').astype(int)

## select useful columns: all columns except not useful ones
df = df[df.columns.difference(['ID', 'INFO', 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT'])]

# Get alleles for each sample
def get_alleles(df):
    for i in df.columns.difference(['CHROM', 'POS']):
        suffix=  str(i) + '_genotype'
        df[suffix] = df[str(i)].astype(str).str[0:3]

get_alleles(df)

## remove original genotype columns
filter_col = [col for col in df if col.endswith('genotype')]
filter_col.append('CHROM')
filter_col.append('POS')
df = df[filter_col]

# replace genotypes: 1/1 by 1, else by 0
list_values = ['0/0', './.', './0', '0/.', '1/0', '0/1']
df = df.replace(to_replace =list_values, value ='NaN')
df = df.replace(to_replace ='1/1', value =1)

a = "Chromosome"
chrom=a + str(df['CHROM'].unique()[0])

## keep sample columns only
df = df[df.columns.difference(['CHROM'])]

# plot SNP density for each sample ==========================================
df = df.set_index('POS').rolling(window_size).count().reset_index().iloc[::increment_value, :]
df = df.melt(id_vars='POS', value_vars=df.columns.difference(['POS']).to_list(), value_name='polym', var_name='sample')

plt.figure(figsize=(9,5))
sns_plot = sns.lineplot(data=df, x='POS',y='polym',hue='sample')
sns_plot.set(xlabel='Chromosome position', ylabel='Number of polymorphisms')
sns_plot.figure.savefig("output.png", dpi=199)
#plt.show()
