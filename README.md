# Spreading-Correction

Supplementary information to the submitted manuscript  "Computational correction of cross-contamination due to index switching" by Anton JM Larsson, Geoff	Stanley, Rahul	Sinha, Irving	L.	Weissman,	Rickard	Sandberg

The Jupyter notebook (Correcting_Spreading_of_Signal_Notebook.ipynb) contains Python code for the analysis and correction of index-swapping, including the generation of Figures 1, 2A-C, S1 and S2 in the manuscript (written by Anton JM Larsson). The notebook (sandbergCorrection_analyzeClustering.ipynb) contains the R code to reproduce Figures 2D-E and S3 of the manuscript (written by Geoff Stanley).

## unspread.py usage example

usage: unspread.py [-h] [--i5 STRING] [--i7 STRING] [--rows INTEGER]
                   [--cols INTEGER] [--idx_col INTEGER] [--sep CHAR]
                   [--h INTEGER] [--c INTEGER] [--t FLOAT]
                   [--idx_in_id BOOLEAN] [--delim_idx CHAR] [--column BOOLEAN]
                   filename

Unspread: Computational correction of barcode index spreading

**positional arguments:**

  filename             .csv file with counts

**optional arguments:**

  -h, --help           show this help message and exit
  
  --i5 STRING          Index name of i5 barcodes (default: 'i5.index.name')
  
  --i7 STRING          Index name of i7 barcodes (default: 'i7.index.name')
  
  --rows INTEGER       Number of rows in plate (default: 16)
  
  --cols INTEGER       Number of columns in plate (default: 24)
  
  --idx_col INTEGER    Which column serves as the index (default: 0)
  
  --sep CHAR           The separator in the .csv file (default ',')
  
  --h INTEGER          The number of reads to use to be considered highly
                       expressed in only one cell (default: 30)
                       
  --c INTEGER          Cutoff to remove addition false positives (default: 5)
  
  --t FLOAT            Threshold for acceptable fraction of spread counts
                       (default: 0.05)
                       
  --idx_in_id BOOLEAN  If the index is in the cell id (i.e. cellid_i5_i7)
                       (Default: 0 (False), set to 1 otherwise (True))
                       
  --delim_idx CHAR     If the index is in the cell id, the delimiting
                       character (Default: '_')
                       
  --column BOOLEAN     If each column is represents a cell, otherwise each
                       row. (default: 1 (True), set to 0 otherwise (False))

The unspread.py script requires a table of read counts supplied as a .csv file with added information regarding each cell's index barcodes. For example from the first plate in the manuscript:

|cell.name | N.index.name |	S.index.name |	0610005C13Rik |	0610007C21Rik  | ...|
| --- | --- | --- | --- | --- | --- |
|HSC02_a_p1c7r2_P01 |	N701	| S522	| 0	| 117 | ...|
|HSC02_a_p1c5r5_P03 | 	N702	| S522	| 0	| 5	| ...|

In this particular example, genes are structured by column and cells by rows but the converse is also supported.

To run the correction of the first plate in the manuscript:

./unspread.py mHSC_plate1HiSeq_counts_IndexInfo.csv --i5 'S.index.name' --i7 'N.index.name' --column 0 --sep ' '

With expected output:

Reading file: mHSC_plate1HiSeq_counts_IndexInfo.csv

Estimating spreading from mHSC_plate1HiSeq_counts_IndexInfo.csv

Found expression to be biased along a certain column and row combination 753 times out of 899

Estimated the median rate of spreading to be 0.0098

Estimated fraction of spread reads to be 0.14827 and variance explained R-squared = 0.8996

Saving figure from analysis to mHSC_plate1HiSeq_counts_IndexInfo_figures.pdf

Saving log file from analysis to mHSC_plate1HiSeq_counts_IndexInfo_unspread.log

Correcting spreading for each gene

Saving correction to mHSC_plate1HiSeq_counts_IndexInfo_corrected.csv
