# 08.18.20
# check coverage at variant positions
# Usage: check_coverage.py <dir with results files> <coverage> <percent positions>

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

# colnames
# CHROM	POS	REF	ALT	$sample.GT	$sample.DP	$sample.AD

# assign arguments
input_dir, coverage, num_positions = sys.argv[1:4]
# initialize counter
index_counter = 0
# initialize results df
results_df = pd.DataFrame()
# loop over files
files = [f for f in os.listdir(input_dir) if f.endswith('annotated.tsv')]
files_sort = sorted(files)
for x in files_sort:
	# get sample name from file name
	name = x.replace('.annotated.tsv', '')
	# open pd df
	df = pd.read_csv(x, sep='\t')
	# filter df by coverage
	cov = df.loc[df[f'{name}.DP'] >= int(coverage)]
	# calculate percent coverage
	cov_percent = round(len(cov) / len(df) * 100, 2)
	# add to results df
	cov_df = pd.DataFrame.from_dict({index_counter: {"sample": name, "coverage_percent": cov_percent}}, orient='index')
	results_df = results_df.append(cov_df)
	index_counter += 1

# write final df to file
results_df.to_csv(f"sample_coverage_{coverage}.csv", index=False)

# how many sample have sufficient coverage
percent_df = results_df.loc[results_df['coverage_percent'] > int(num_positions)]
sample_percent = round(len(percent_df) / len(files_sort) * 100, 1)
# print results
print(f'{sample_percent}% of samples achieved >= {coverage}x at >{num_positions}% of positions captured by the custom reagent')

# create histogram of results
plt.hist(results_df['coverage_percent'], bins=25)
plt.ylabel("Sample count")
plt.xlabel(f'% positions with {coverage} coverage')
plt.show()
plt.savefig(f"cov_histogram_{coverage}.png", format='png')
	
