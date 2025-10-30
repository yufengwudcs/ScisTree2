# Install for Python

## Required Tools

To use ScisTree2, you will need the following tools and libraries installed:
*   **`python` & `pip`**: Version 3.6 or higher.
*   **`g++`**: A C++ compiler. (from source only)
*   **`make`**: For building the backend. (from source only)

```{note}
We have successfully tested it on Linux, macOS, and Windows (via WSL).
```

## Installation
### 1. Install from [PyPI](https://pypi.org/project/scistree2/):
1. Upgrade pip:
    ```bash
    python -m pip install --upgrade pip
    ```
2. Install `scistree2`:
    ```bash 
    pip install scistree2
    ```


### 2. Install from source:
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

    ```{note}
    We recommend that users create a virtual environment using either `conda` or `venv` to comply with [PEP 668](https://peps.python.org/pep-0668/).
    ```
    <!-- Alternatively, you can use `python setup.py install`:
    ```bash
    python setup.py install
    ```
    The `setup.py` script is configured to first build the C++ executable (similar to running `make` in the `src` directory) and then include it in the Python package. -->
