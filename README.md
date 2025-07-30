<p align="center">
  <img src="img/logo.png" alt="ScisTree2 Logo" width="200">
</p>

<h2 align="center">Fast cell lineage tree reconstruction and genotype calling for large single cell DNA sequencing data.</h2>

<!-- ## Introduction -->
Software accompanyment for [*Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single-Cell Data Using Efficient Local Search*, Haotian Zhang, Yiming Zhang, Teng Gao and Yufeng Wu, manuscript, 2025](https://www.biorxiv.org/content/10.1101/2024.11.08.622704v1) (under the title *"ScisTree2: An Improved Method for Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single Cell Data"*). This work was presented in the RECOMB 2025 conference. The ScisTree2 paper is currently under review.

This is an enhanced version of ScisTree (*Accurate and efficient cell lineage tree inference from noisy single cell data: the maximum likelihood perfect phylogeny approach, Bioinformatics, Wu, Volume 36, Issue 3, Pages 742–750, 2020*).

## Required Tools

To use ScisTree2, you will need the following tools and libraries installed:
*   **`python` & `pip`**: Version 3.6 or higher.
*   **`g++`**: A C++ compiler.
*   **`make`**: For building the backend.

**We have successfully tested it on Linux, macOS, and Windows (via WSL).*

## Installation/

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/yufengwudcs/ScisTree2.git
    cd scistree2
    ```

2.  **Install the Python package (includes C++ backend compilation):**
    You can install the `scistree2` package using `pip`:
    ```bash
    pip install .
    ```
    This command will also automatically compile the C++ backend. Once built, the executable binary file can be found in `scistree2/bin`.

    **We recommend that users create a virtual environment using either `conda` or `venv` to comply with [PEP 668](https://peps.python.org/pep-0668/).*
    <!-- Alternatively, you can use `python setup.py install`:
    ```bash
    python setup.py install
    ```
    The `setup.py` script is configured to first build the C++ executable (similar to running `make` in the `src` directory) and then include it in the Python package. -->

3.  **(Optional) Manual C++ backend build (for testing/development):**
    If you want to build or test the C++ backend (`scistree`) independently, you can navigate to the `src` directory and compile it using `make`:
    ```bash
    cd src
    make
    # You can then test it directly, e.g., ./scistree triv4-paper-1.txt
    ```
    This step is not required for the Python package installation if using `pip install .` as described above. See more details below.

## Tutorial

A detailed tutorial on how to use ScisTree2 is available as a Jupyter Notebook in the `tutorials/` directory:

*   **[ScisTree2 Tutorial](tutorials/ScisTree2_Tutorial.ipynb)**

The tutorial covers:
*   Getting started with ScisTree2.
*   Running inference with probabilistic genotype matrices.
*   Running inference with raw read data.
*   Running inference with VCF file.
*   Visualizing trees.
*   Evaluating results using various metrics.

The example data used in the tutorial can be found in the `tutorials/data/` directory.

## Using ScisTree2 by mannually building the C++ code.
Download the source code from GitHub repository. Decompress it if you download as a zip file. Open a console window and enter the main source code directory called "src". Type "make". That should be all you need!

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

### Data format of ScisTree2?
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

* -T <num-of-threads>:  ScisTree2 now supports multi-threading.
* -s <max-num-of-iterations>: To control the running time, you can specifiy the maximum number of iterations to run. Specify a smaller number (e.g., 5) can reduce the running time (but may reduce the accuracy). Default: 1,000 iterations.

#### You may also read the ScisTree2's User Manual, which is in PDF format and is distributed as part of ScisTree2. 

## What is new about ScisTree2 over ScisTree1?
The main change is about speed and accuracy. ScisTree2 is order of mangnitude faster than ScisTree1. ScisTree2 supports multi-threading while ScisTree1 doesn't. More importantly, ScisTree2 implements faster and also possibly more accurate tree search algorithms. By default, ScisTree2 performs the subtree prune and regraft (SPR) local search, while ScisTree1 performs neareast neighbor interchange (NNI) search. The SPR local search is usually more accurate than the NNI search. Our tests show that ScisTree2 can infer cell lineage tree from data with 10,000 cells (and say 10,000 single nucleiotide variant or SNV sites) while being more accurate in both cell lineage tree and genotype calling. 

## Data Availability 

All simulated data, experimental data(HGSOC), and scripts used to reproduce the results in paper *"Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single-Cell Data Using Efficient Local Search"* are released at Zenodo.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15620911.svg)](https://zenodo.org/records/15620911)

## Contact
Post your issues here inside GitHub repositary if you have questions/issues.
