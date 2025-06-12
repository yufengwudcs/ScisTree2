<p align="center">
  <img src="img/logo.png" alt="SciSTree2 Logo" width="200">
</p>

<h2 align="center">Fast cell lineage tree reconstruction and genotype calling for large single cell DNA sequencing data.</h2>

<!-- ## Introduction -->
Software accompanyment for [*Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single-Cell Data Using Efficient Local Search*, Haotian Zhang, Yiming Zhang, Teng Gao and Yufeng Wu, manuscript, 2025](https://www.biorxiv.org/content/10.1101/2024.11.08.622704v1) (under the title *"ScisTree2: An Improved Method for Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single Cell Data"*). This work was presented in the RECOMB 2025 conference. The ScisTree2 paper is currently under review.

This is an enhanced version of ScisTree (*Accurate and efficient cell lineage tree inference from noisy single cell data: the maximum likelihood perfect phylogeny approach, Bioinformatics, Wu, Volume 36, Issue 3, Pages 742–750, 2020*).

## Required Tools

To use SciSTree2, you will need the following tools and libraries installed:

*   **python & pip**: Version 3.6 or higher.
*   **g++**: A C++ compiler (for building the backend).
*   **make**: The `make` utility (for building the backend).

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/yufengwudcs/ScisTree2.git
    cd ScisTree2
    ```

2.  **Install the Python package (includes C++ backend compilation):**
    You can install the `scistree2` package using `pip`:
    ```bash
    pip install .
    ```
    This command will also automatically compile the C++ backend. Once built, the executable binary file can be found in `scistree2/bin`.

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
    This step is not required for the Python package installation if using `pip install .` as described above.

## Tutorial

A detailed tutorial on how to use SciSTree2 is available as a Jupyter Notebook in the `tutorials/` directory:

*   **[SciSTree2 Tutorial](tutorials/Scistree2_Tutorial.ipynb)**

The tutorial covers:
*   Getting started with SciSTree2.
*   Running inference with probabilistic genotype matrices.
*   Running inference with raw read data.
*   Visualizing trees.
*   Evaluating results using various metrics.

The example data used in the tutorial can be found in the `tutorials/data/` directory.

## Data Availability 

All simulated data, experimental data(HGSOC), and scripts used to reproduce the results in paper *"Large-scale Inference of Cell Lineage Trees and Genotype Calling from Noisy Single-Cell Data Using Efficient Local Search"* are released at Zenodo.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15620911.svg)](https://zenodo.org/records/15620911)

## Contact
Post your issues here inside GitHub repositary if you have questions/issues.
