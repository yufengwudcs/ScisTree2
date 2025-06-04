# ScisTree2
Fast cell lineage tree reconstruction and genotype calling for large single cell DNA sequencing data.  

Current version: v2.2.0.0. Released: October 24, 2024.

Software accompanyment for "Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single-Cell Data Using Efficient Local Search", Haotian Zhang, Yiming Zhang, Teng Gao and Yufeng Wu, manuscript, 2025. The preprint of this paper is at: https://www.biorxiv.org/content/10.1101/2024.11.08.622704v1 (under the title "ScisTree2: An Improved Method for Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single Cell Data"). This work was presented in the RECOMB 2025 conference. The ScisTree2 paper is currently under review.

This is an enhanced version of Scistree1 (Wu, Accurate and efficient cell lineage tree inference from noisy single cell data: the maximum likelihood perfect phylogeny approach, Bioinformatics, 2020). 

# Getting started
Download the source code from GitHub repository. Decompress it if you download as a zip file. Open a console window and enter the main source code directory called "ScisTree2-source-code". Type "make". That should be all you need!

The executable is called "scistree". You can find whether it is built by doing a "ls". 

Check if ScisTree2 is ready to run by typing: "./scistree". You should see some output about the basic usage of ScisTree2. 

Now type: "./scistree triv4-paper-1.txt"
You should see the following output:

*** SCISTREE ver. 2.2.0.0, October 24, 2024 ***   

Called genotypes output to file: Test/triv4-paper-1.txt.genos.imp

**** Maximum log-likelihood: -6.27126, number of changed genotypes: 2

Computed log-lielihood from changed genotypes: -6.27126

Constructed single cell phylogeny: (((1,3),(2,4)),5)

Elapsed time = 0 seconds.

# What is new about ScisTree2 over ScisTree1?
The main change is about speed and accuracy. ScisTree2 is order of mangnitude faster than ScisTree1. ScisTree2 supports multi-threading while ScisTree1 doesn't. More importantly, ScisTree2 implements faster and also possibly more accurate tree search algorithms. By default, ScisTree2 performs the subtree prune and regraft (SPR) local search, while ScisTree1 performs neareast neighbor interchange (NNI) search. The SPR local search is usually more accurate than the NNI search. Our tests show that ScisTree2 can infer cell lineage tree from data with 10,000 cells (and say 10,000 single nucleiotide variant or SNV sites) while being more accurate in both cell lineage tree and genotype calling. 

# How to use ScisTree2?
First, you should understand some basics about ScisTree2. I would recommend to read the user mannual of the orogianl ScisTree: https://github.com/yufengwudcs/ScisTree/blob/master/ScisTree-UserManual.pdf

The first thing to use ScisTree2 is to prepare the input. ScisTree2 uses the same data format as ScisTree1. Here is the content of triv4-paper-1.dat:

HAPLOID 6 5  

0.01 0.6 0.08 0.8 0.7   

0.8 0.02 0.7 0.01 0.3   

0.02 0.8 0.02 0.8 0.9   

0.9 0.9 0.8 0.8 0.02   

0.01 0.8 0.01 0.8 0.9   

0.05 0.02 0.7 0.05 0.9  


* Explanations. HAPLOID: specify binary input (at the moment this is the only format supported). 6: number of SNV sites; 5: number of cells; each following row: the probability of the five cells being zero (wild-type). **Be careful: the rows are for the SNV sites and the columns are for the cells. Don't get this wrong.**

ScisTree2 is essentially a faster and also somewhat more accurate ScisTree. Some features from the original ScisTree (version 1) are not supported in the current implementaiton of ScisTree2. These include: (i) ternary data input: ScisTree2 only supports binary data as of now; (ii) parameter imputation and doublet imputation. I haven't got chance to upgrade these features. For the moment, ScisTree2 is dedicated for cell lineage tree inference.

The following options can be useful.

	 -e                Output mutation tree (may not be binary tree) from called genotypes branch labels.
  
	 -e0               Output mutation tree but don't output labels (for visualizing large trees).
  
	 -q                Use NNI local tree search (NNI is faster but less accurate). By default, ScisTree2 uses SPR local search. Our experience shows the default SPR search is usually very fast. Still, you can try the NNI local search if you want.

There are options that are new to ScisTree2.

* -T num-of-threads:  ScisTree2 now supports multi-threading. 

## You may also read the ScisTree2's User Manual, which is in PDF format and is distributed as part of ScisTree2. 

## To learn more about ScisTree2, see our [Scistree2 Tutorial](https://github.com/yufengwudcs/ScisTree2/blob/main/Scistree2_Tutorial.ipynb).

# Paper Data Availability
All simulated data, real data(HGSOC), and scirpts used to reproduce the results in paper "ScisTree2: An Improved Method for Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single Cell Data. 10.1101/2024.11.08.622704. Zhang, H & Zhang, Y & Gao, T & Wu, Y. (2024)." can be found in Github Release.

# Contact
Post your issues here inside GitHub repositary if you have questions/issues.
