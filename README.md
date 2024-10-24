# ScisTree2
Fast cell lineage tree reconstruction for large single cell data.  

Current version: v2.1.0.0. Released: January 25, 2024.

Software accompanyment for "ScisTree2: An Improved Method for Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single Cell Data", Haotian Zhang, Yiming Zhang, Teng Gao and Yufeng Wu, submitted for publication, 2024.

This is an enhanced version of Scistree1 (Wu, Accurate and efficient cell lineage tree inference from noisy single cell data: the maximum likelihood perfect phylogeny approach, Bioinformatics, 2020). 

# Getting started
Download the source code (the zip file in this repositary). Decompress it. Open a console window and go the main source code directory. Type "make". That should be all you need!

Check if ScisTree2 is ready to run by typing: "./scistree". You should see some output about the basic usage of ScisTree2. 

Now type: "./scistree example-input.dat"
You should see the following output:

*** SCISTREE ver. 2.1.0.0, Janurary 25, 2024 ***

**** Maximum log-likelihood: -6.27126, number of changed genotypes: 2
Computed log-lielihood from changed genotypes: -6.27126
Constructed single cell phylogeny: (((1,3),(2,4)),5)
Elapsed time = 0 seconds.

# What is new about ScisTree2 over ScisTree1?
The main change is about speed. ScisTree2 is order of mangnitude faster than ScisTree1. ScisTree2 supports multi-threading while ScisTree1 doesn't. More importantly, ScisTree2 implements faster and also possibly more accurate tree search algorithms. Our tests show that ScisTree2 can infer cell lineage tree from data with 10,000 cells (and say 10,000 single nucleiotide variant or SNV sites). 

# How to use ScisTree2?
First, you should understand some basics about ScisTree2. I would recommend to read the user mannual of the orogianl ScisTree: https://github.com/yufengwudcs/ScisTree/blob/master/ScisTree-UserManual.pdf

The first thing to use ScisTree2 is to prepare the input. ScisTree2 uses the same data format as ScisTree1. Here is the content of example-input.dat:

HAPLOID 6 5
0.01 0.6 0.08 0.8 0.7
0.8 0.02 0.7 0.01 0.3
0.02 0.8 0.02 0.8 0.9
0.9 0.9 0.8 0.8 0.02
0.01 0.8 0.01 0.8 0.9
0.05 0.02 0.7 0.05 0.9

* Explanations. HAPLOID: specify binary input (at the moment this is the only format supported). 6: number of SNV sites; 5: number of cells; each following row: the probability of the five cells being zero (wild-type).

ScisTree2 is essentially a faster and also somewhat more accurate ScisTree. Some features from the original ScisTree (version 1) are not supported in the current implementaiton of ScisTree2. These include: (i) ternary data input: ScisTree2 only supports binary data as of now; (ii) parameter imputation and doublet imputation. I haven't got chance to upgrade these features. For the moment, ScisTree2 is dedicated for cell lineage tree inference.

There are options that are new to ScisTree2.

* -s: turn on SPR local search. Example: "./scistree -s example-input.dat". This performs SPR local search which is somewhat slower than the default NNI search, but is usually more accurate. Our tests show that the SPR mode can run for data with as many as 10,000 cells.

# Contact
Post your issues here inside GitHub repositary if you have questions/issues.
