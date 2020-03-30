This repository contains scripts and description of methods for the  RNA sequencing analyses for the manuscript Combining experimental evolution and genomics to understand how seed beetles adapt to a marginal host plant (Authors: Alex Rego, Samridhi Chaturvedi, Amy Springer, Alexandra Lish, Caroline Barton, Karen Kapheim, Frank Messina, Zachariah Gompert). March 2020

Here are the details of the files in the repository:

1. sequence_alignments.md = Pipeline to do quality filtering, sequence alignment, and creating counts matrix for RNA sequencing data.
2. 
3. differential_expression_analyses.R = R code to do count normalization and differential expression analyses.
4. prepare_countmatrix.R = R code to prepare a count matrix by combining data from all samples in the analyses.
5. transcript_annot.py = Python script to combine count matrix with gene IDs from genome annotation (gff) file.
