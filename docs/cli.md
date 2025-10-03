# Using ScisTree2 in C++

To run ScisTree2 directly from the console, please refer to [C++ installation](./installation/installation-c.md) in the previous chapter.

The executable is called `scistree`. Check if ScisTree2 is ready to run by typing: `./scistree`, you should see some output about the basic usage of ScisTree2.

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

```{note}
Be careful: the rows are for the SNV sites and the columns are for the cells. Don't get this wrong.
```