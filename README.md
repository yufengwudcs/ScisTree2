<p align="center">
  <img src="img/logo.png" alt="ScisTree2 Logo" width="200">
</p>

<h2 align="center">Fast cell lineage tree reconstruction and genotype calling for large single cell DNA sequencing data </h2>

<!-- ## Introduction -->
Software accompanyment for **"ScisTree2 enables large-scale inference of cell lineage trees and genotype calling using efficient local search"**, Haotian Zhang, Yiming Zhang, Teng Gao and Yufeng Wu, Genome Research, in press, 2025 (https://genome.cshlp.org/content/early/2025/09/03/gr.280542.125.abstract). Here is the preprint in a slightly different title [*Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single-Cell Data Using Efficient Local Search*, Haotian Zhang, Yiming Zhang, Teng Gao and Yufeng Wu, manuscript, 2025](https://www.biorxiv.org/content/10.1101/2024.11.08.622704v1) (under the title *"ScisTree2: An Improved Method for Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single Cell Data"*). This work was presented in the RECOMB 2025 conference.

This is an enhanced version of ScisTree (*Accurate and efficient cell lineage tree inference from noisy single cell data: the maximum likelihood perfect phylogeny approach, Bioinformatics, Wu, Volume 36, Issue 3, Pages 742–750, 2020*).

## Required Tools

To use ScisTree2, you will need the following tools and libraries installed:
*   **`python` & `pip`**: Version 3.6 or higher.
*   **`g++`**: A C++ compiler.
*   **`make`**: For building the backend.

**We have successfully tested it on Linux, macOS, and Windows (via WSL).*

## Installation

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
    # You can then test it directly, e.g., ./scistree example_input.txt
    ```
    This step is not required for the Python package installation if using `pip install .` as described above. See more details below.

## Tutorial
ScisTree2 offers both Python and C++ interfaces. We recommend using the Python version because it provides a wider variety of supported input formats and evaluation tools, and it is more easily integrated into the broader Python ecosystem.

### Using ScisTree2 in Python
A detailed tutorial on how to use ScisTree2 in Python is available as a Jupyter Notebook in the `tutorials/` directory:

*   **[ScisTree2 Tutorial](tutorials/ScisTree2_Tutorial.ipynb)**

The tutorial covers:
*   Getting started with ScisTree2.
*   Running inference with probabilistic genotype matrices (CSV supported).
*   Running inference with raw read data (CSV supported).
*   Running inference with VCF file.
*   Visualizing trees.
*   Evaluating results using various metrics.
*   Bootstrapping for branch (clade) confidence estimates (added September 27, 2025). 

The example data used in the tutorial can be found in the `tutorials/data/` directory.

### Using ScisTree2 in C++
To run ScisTree2 directly from the console, please refer to step 3 in the [installation](#installation) guide above.

The executable is called `scistree`. 
Check if ScisTree2 is ready to run by typing: `./scistree`, you should see some output about the basic usage of ScisTree2. 

Now type: 
`./scistree example_input.txt`, you should see the following output:
```text
*** SCISTREE ver. 2.2.3.0, August 14, 2025 ***

#cells: 5, #sites: 6
List of cell names: c1 c2 c3 c4 c5 
Called genotypes output to file: example_input.txt.genos.imp
**** Maximum log-likelihood: -6.27126, number of changed genotypes: 2
Computed log-lielihood from changed genotypes: -6.27126
Constructed single cell phylogeny: (((c1,c3),(c2,c4)),c5)
Elapsed time = 0 seconds.
```

**Options:**
* `-e`: Output a mutation tree (which may not be a binary tree) with branch labels from the called genotypes.
* `-e0`: Output a mutation tree without branch labels, which is useful for visualizing large trees.
* `-q`: Use NNI (Nearest Neighbor Interchange) for local tree search. NNI is faster but less accurate. By default, ScisTree2 uses SPR (Subtree Pruning and Regrafting) local search, which we have found to be very fast.
* `-T <num-of-threads>`: Specify the number of threads for multi-threading support.
* `-s <num-of-iterations>`: Set the maximum number of iterations to control the running time. A smaller number (e.g., 5) will reduce the running time but may also reduce accuracy. Default: 1,000 iterations.

You may also read the ScisTree2's User Manual, which is in PDF format and is distributed as part of ScisTree2. 
#### Data format of ScisTree2 in C++?
First, you should understand some basics about ScisTree2. I would recommend to read the [user mannual](https://github.com/yufengwudcs/ScisTree/blob/master/ScisTree-UserManual.pdf) of the orgianl ScisTree.

The first thing to use ScisTree2 is to prepare the input. Here is the content of an example(example_input.txt):
```js
c1 c2 c3 c4 c5
s1 0.01 0.6 0.08 0.8 0.7
s2 0.8 0.02 0.7 0.01 0.3
s3 0.02 0.8 0.02 0.8 0.9
s4 0.9 0.9 0.8 0.8 0.02
s5 0.01 0.8 0.01 0.8 0.9
s6 0.05 0.02 0.7 0.05 0.9
```

Explanations: 
- You should specifiy the cell names in the first row. For example, "c1 c2 c3 c4 c5". Please note that don't use **HAPLOID** or **HAPLOTYPES** as cell names. These two words are reserved keywords in ScisTree2.
- The following row starts with the row identifier, then the probability of the five cells being zero (wild-type). For example, the second row says for the first site, the probability of the first cell (cell 1) has probability 0.01 being the wild type, the second cell has probability 0.6 being the wild type, and so on.

    **Be careful: the rows are for the SNV sites and the columns are for the cells. Don't get this wrong.**

ScisTree2 is essentially a faster and also somewhat more accurate ScisTree. Some features from the original ScisTree (version 1) are not supported in the current implementaiton of ScisTree2. These include: (i) ternary data input: ScisTree2 only supports binary data as of now; (ii) parameter imputation and doublet imputation. I haven't got chance to upgrade these features. For the moment, ScisTree2 is dedicated for cell lineage tree inference.

## What is new about ScisTree2 over ScisTree?
The main change is about speed and accuracy. ScisTree2 is order of mangnitude faster than ScisTree. ScisTree2 supports multi-threading while ScisTree doesn't. More importantly, ScisTree2 implements faster and also possibly more accurate tree search algorithms. By default, ScisTree2 performs the subtree prune and regraft (SPR) local search, while ScisTree performs neareast neighbor interchange (NNI) search. The SPR local search is usually more accurate than the NNI search. Our tests show that ScisTree2 can infer cell lineage tree from data with 10,000 cells (and say 10,000 single nucleiotide variant or SNV sites) while being more accurate in both cell lineage tree and genotype calling. 

## Data Availability 

All simulated data, experimental data(HGSOC), and scripts used to reproduce the results in the SicsTree2 paper are released at Zenodo.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15620911.svg)](https://zenodo.org/records/15620911)

## Contact
Post your issues here inside GitHub repositary if you have questions/issues.
