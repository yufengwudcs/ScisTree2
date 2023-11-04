# ScisTree2
Fast cell lineage tree reconstruction

Software accompanyment for "ScisTree2: Accurate Cell Lineage Tree Inference for Tens of Thousands of Cells from Noisy Single Cell Data", Haotian Zhang, Teng Gao, Yiming Zhang, Peter Kharchenko and Yufeng Wu, submitted for publication, 2023.

This is an enhanced version of Scistree (Wu, Accurate and efficient cell lineage tree inference from noisy single cell data: the maximum likelihood perfect phylogeny approach, Bioinformatics, 2020). 

# Getting started
Download the source code (the zip file in this repositary). Decompress it. Open a console window and go the main source code directory. Type "make". That should be all you need!

Check if ScisTree2 is ready to run by typing: "./scistree". You should see some output about the basic usage of ScisTree2. 

Now type: "./scistree example-input.dat"
You should see the following output:

*** SCISTREE ver. 2.0.0.0, November 4, 2023 ***

**** Maximum log-likelihood: -6.27126, number of changed genotypes: 2
Computed log-lielihood from changed genotypes: -6.27126
Constructed single cell phylogeny: (((1,3),(2,4)),5)
Elapsed time = 0 seconds.

# How to use ScisTree2?
First, you should understand some basics about ScisTree2. 
