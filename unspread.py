#!/usr/bin/env python3
import os, sys,re
import argparse
import pandas as pd
import numpy as np
from matplotlib import use
import matplotlib.pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests
import statsmodels.api as sm
from scipy.stats import hypergeom
from scipy.linalg import solve_sylvester
from scipy.stats import binom_test
from patsy import dmatrices

# Used for getting genes to estimate spreading
def num_counts(col, high):
    return np.sum(col > high)


parser = argparse.ArgumentParser(description='Unspread: Computational correction of barcode index spreading')
parser.add_argument('filename', metavar='filename', type=str, nargs=1,help='.csv file with counts' )
parser.add_argument('--i5', metavar='STRING', type=str, nargs=1, default=['i5.index.name'], help='Index name of i5 barcodes (default: \'i5.index.name\')')
parser.add_argument('--i7', metavar='STRING', type=str, nargs=1, default=['i7.index.name'], help='Index name of i7 barcodes  (default: \'i7.index.name\')')
parser.add_argument('--rows', metavar='INTEGER', default=[16], type=int, nargs=1, help='Number of rows in plate (default: 16)')
parser.add_argument('--cols', metavar='INTEGER', default=[24], type=int, nargs=1,  help='Number of columns in plate (default: 24)')
parser.add_argument('--idx_col', metavar='INTEGER', default=[0], type=int, nargs=1,  help='Which column serves as the index (default: 0)')
parser.add_argument('--sep', metavar='CHAR', default=[','], type=str, nargs=1,  help='The separator in the .csv file (default \',\')')
parser.add_argument('--h', metavar='INTEGER', default=[30], type=int, nargs=1,  help='The number of reads to use to be considered highly expressed in only one cell (default: 30)')
parser.add_argument('--c', metavar='INTEGER', default=[5], type=int, nargs=1, help='Cutoff to remove addition false positives (default: 5)')
parser.add_argument('--t', metavar='FLOAT', default=[0.05], type=float, nargs=1, help='Threshold for acceptable fraction of spread counts (default: 0.05)')
parser.add_argument('--idx_in_id', metavar='BOOLEAN', default=[0], type=float, nargs=1, help='If the index is in the cell id (i.e. cellid_i5_i7) (Default: 0 (False), set to 1 otherwise (True))')
parser.add_argument('--delim_idx', metavar='CHAR', default=['_'], type=str, nargs=1, help='If the index is in the cell id, the delimiting character (Default: \'_\')')
parser.add_argument('--column', metavar='BOOLEAN', default=[1], type=int, nargs=1, help='If each column is represents a cell, otherwise each row. (default: 1 (True), set to 0 otherwise (False))')
parser.add_argument('--rate', metavar='FLOAT', default=[0.0], type=float, nargs=1, help='Set spreading rate manually (overrides any estimated rate).')

args = parser.parse_args()
filename = args.filename[0]
i5_index_name = args.i5[0]
i7_index_name = args.i7[0]
separator = args.sep[0]
n_rows = args.rows[0]
n_cols = args.cols[0]
idx = args.idx_col[0]
high = args.h[0]
c = args.c[0]
threshold = args.t[0]
column = args.column[0]
idx_in_id = args.idx_in_id[0]
delim_idx = args.delim_idx[0]
r = args.rate[0]

print('Reading file: {}'.format(filename))

# Load in count table
df = pd.read_csv(filename, index_col=idx, sep=separator)
# Transform the data frame if each column represents a cell
if np.bool_(column):
    df = df.T
# Turns the end of the cell id strings into index name ids
if np.bool_(idx_in_id):
    i5_index_list = []
    i7_index_list= []
    for i,string in enumerate(df.index.values):
        sp = string.split(delim_idx)
        i5_index_list.append(sp[-2])
        i7_index_list.append(sp[-1])
    df[i5_index_name] = i5_index_list
    df[i7_index_name] = i7_index_list
df = df.sort_values(by=[i5_index_name,i7_index_name], ascending=[False, True])
df = df.loc[:,~df.columns.duplicated()]
df_noindex = df.drop([i5_index_name,i7_index_name], axis=1).astype(np.int)

# Check if the shape makes sense
if (df_noindex.shape[0] != n_rows*n_cols):
    print('Number of cells in count file not the same as specified, exiting...')
    quit()

base = os.path.splitext(os.path.basename(filename))[0]

print('Estimating spreading from {}'.format(filename))

# Get genes to estimate the rate of spreading and fraction of contaminating reads.
names_list = []
for i, col in df_noindex.items():
    if num_counts(col, high) == 1:
        names_list.append(i)

n_names = len(names_list)
if n_names == 0:
    print('Found no genes useable to estimate spreading, exiting...')
    quit()
# Count the number of true and spread counts respectively
true_counts = np.zeros(n_names)
spread_counts = np.zeros((n_names, n_rows+n_cols - 2))
num_wells_counts = np.zeros(n_names)

for i,namn in enumerate(names_list):
    if len(df_noindex[namn].values) == n_rows*n_cols:
        true_counts[i] = np.amax(df_noindex[namn].values)
        w = np.argwhere(df_noindex[namn].values.reshape(n_rows,n_cols) == true_counts[i])[0]
        spread_counts[i] = np.append(np.append(df_noindex[namn].values.reshape(n_rows,n_cols)[w[0], :w[1]],df_noindex[namn].values.reshape(n_rows,n_cols)[w[0], w[1]+1:] ), np.append(df_noindex[namn].values.reshape(n_rows,n_cols)[:w[0], w[1]],df_noindex[namn].values.reshape(n_rows,n_cols)[w[0]+1:, w[1]] ))
        num_wells_counts[i] = np.sum(df_noindex[namn].values != 0) - 1


prop_spread = np.zeros((len(names_list), n_rows+n_cols - 2))

for i in range(len(names_list)):
    prop_spread[i] = spread_counts[i]/true_counts[i]

test_bias = np.zeros(n_names)

# Test whether the spread counts are biased along the column and row of the source well.
for i in range(len(names_list)):
    test_bias[i] = hypergeom.sf(np.sum(spread_counts[i] != 0)-1, n_rows*n_cols-1, n_rows+n_cols - 2, num_wells_counts[i])

mt = multipletests(test_bias[test_bias != 1])[0]
bias_log = 'Found expression to be biased along a certain column and row combination {} times out of {}'.format(np.sum(mt), len(mt))
print(bias_log)
f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7,14))
ax1.hist(test_bias)
ax1.set_title('Bias Test')
ax1.set_xlabel('p value counts bias')
ax1.set_ylabel('Frequency (genes)')

if np.sum(mt) != 0:

    rate_spreading = np.median(prop_spread[test_bias != 1][mt].flatten()[prop_spread[test_bias != 1][mt].flatten() != 0])
    rate_log = 'Estimated the median rate of spreading to be {}'.format(np.round(rate_spreading,5))
    print(rate_log)

    ax2.hist(np.log10(prop_spread[test_bias != 1][mt].flatten()[prop_spread[test_bias != 1][mt].flatten() != 0]))
    ax2.set_title('Rate of Spreading')
    ax2.text(0.7, 0.9,'Median: ' + str(np.round(rate_spreading,5)), ha='center', va='center', transform=ax2.transAxes)
    ax2.set_ylabel('Frequency (Wells)')
    ax2.set_xlabel(r'Rate of spreading ($log_{10}$)')

    # Use linear regression to estimate the fraction of spread reads
    spread_counts = spread_counts[test_bias != 1][mt]
    true_counts = true_counts[test_bias != 1][mt]
    true_median = np.median(true_counts)
    model_df = pd.DataFrame([true_counts, np.apply_along_axis(np.sum, 1,spread_counts)], index=['true', 'spread']).T
    y,X = dmatrices('spread ~  true', data = model_df, return_type='dataframe')
    model = sm.OLS(y, X).fit()
    ols_log = 'Estimated fraction of spread reads to be {} and variance explained R-squared = {}'.format(np.round(model.params['true'],5), np.round(model.rsquared,5))
    print(ols_log)

    ax3.scatter(np.log10(true_counts), np.log10(np.apply_along_axis(np.sum, 1,spread_counts)), s = 1)
    ax3.set_xlabel(r'log$_{10}$(Source counts)')
    ax3.set_ylabel(r'log$_{10}$(Spread counts)')
    ax3.text(0.3, 0.9,r'coef = {}, $R^2 = $ {}'.format(np.round(model.params['true'],5), np.round(model.rsquared,5)), ha='center', va='center', transform=ax3.transAxes)
    ax3.set_title('Genes highly expressed in only one cell')
else:
    print('Found no genes with bias along the column and row combination')
    rate_log = 'Found no genes with bias along the column and row combination'
    ols_log = 'Found no genes with bias along the column and row combination'

print('Saving figure from analysis to {}'.format('{}_figures.pdf'.format(base)))
plt.tight_layout()
plt.savefig('{}_figures.pdf'.format(base))

print('Saving log file from analysis to {}'.format('{}_unspread.log'.format(base)))

if np.sum(mt) != 0:
    with open('{}_unspread.log'.format(base), "w") as log_file:
        print('# {}\n# {}\n# {}\n# {}\n{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(filename, bias_log, rate_log, ols_log, filename, np.sum(mt), len(mt), np.round(rate_spreading,7), np.round(model.params['true'],7), np.round(model.rsquared,7), true_median), file=log_file)
else:
    with open('{}_unspread.log'.format(base), "w") as log_file:
        print('# {}\n# {}\n# {}\n# {}\n{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(filename, bias_log, rate_log, ols_log, filename, np.sum(mt), len(mt), 'NaN', 'NaN', 'NaN', 'Nan'), file=log_file)

# You can set the threshold yourself
if (np.sum(mt) == 0 or model.params['true'] < threshold) and r == 0.0:
    print('The experiment shows no or an acceptable amount of spreading, correction is not neccessary. Exiting...')
    quit()
if r != 0.0:
    rate_spreading = r
# Setting the rate of spreading matrices
column_spread = np.zeros((n_rows, n_rows))
row_spread = np.zeros((n_cols,n_cols))
column_spread[:,:] = rate_spreading
row_spread[:,:] = rate_spreading
np.fill_diagonal(column_spread, 0.5 - rate_spreading/2)
np.fill_diagonal(row_spread, 0.5 - rate_spreading/2)

# The function which does the correction, you can set the cutoff yourself
def adjust_reads(mat, column_spread = column_spread, row_spread = row_spread, cutoff = c, r = n_rows, c = n_cols):
    mat = np.array(mat).flatten()
    mat = mat.reshape(r,c)
    adjusted_reads = np.rint(solve_sylvester(column_spread,row_spread,mat))
    # A lower bound cutoff removes false positives (unfortunately also remove true reads with low counts in that cell)
    adjusted_reads[adjusted_reads < cutoff] = 0
    return adjusted_reads

print('Correcting spreading for each gene')

adj_list = []
for i, col in df_noindex.items():
    adj_list.append(adjust_reads(col).flatten())

# Put correction into a dataframe and save to a .csv file
df_adj = pd.DataFrame(data=adj_list, index = df_noindex.columns.values, columns= df_noindex.index.values).T

df_adj = pd.concat([df[i7_index_name], df[i5_index_name] , df_adj], axis = 1)

print('Saving correction to {}'.format('{}_corrected.csv'.format(base)))

df_adj.to_csv('{}_corrected.csv'.format(base))
