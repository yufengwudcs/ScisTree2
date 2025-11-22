# Here, we provide a Python interface for conveniently calling Scistree2.
import os
import sys
import multiprocessing as mp
import uuid
import numpy as np
import pandas as pd
import subprocess as sp
import importlib.resources
from .treeutils import *


def get_executable_path(provided_path=None):
    if provided_path:
        return provided_path
    binary_name = "scistree.exe" if sys.platform == "win32" else "scistree"
    return os.path.join(os.path.dirname(__file__), 'bin', binary_name)


class ScisTree2:
    """
    Scistree2 Caller.

    Args:
        bin_path: Path to binary executable file.
        threads: Number of threads in use. Default uses all.
        nj: Call NJ only, serves when M is very big.
        spr: Enable SPR local search. Default is on.
        nni: Enable NNI local search. Default is off.
        iterative: Enable iterative optimization. Default is off.
        verbose: Show outputs.
    """

    def __init__(
        self, threads=-1, nj=False, spr=True, nni=False, max_iter=0, verbose=True
    ):
        self.bin_path = get_executable_path()
        self.nj = nj
        self.spr = spr
        self.nni = nni
        assert max_iter >= 0, "max_iter should be positive."
        self.max_iter = max_iter
        self.cmd = self.build_cmd(self.bin_path, threads, nj, nni, verbose)

    def build_cmd(self, bin_path, threads, nj, nni, verbose):
        if threads == -1:
            threads = mp.cpu_count()
        cmd = [bin_path, "-T", str(threads)]
        if verbose:
            cmd.append("-v")
        if nj:
            cmd.append("-n")
        if nni:
            cmd.append("-q")
        if self.max_iter:
            cmd.append("-s")
            cmd.append(str(self.max_iter))
        return cmd

    @staticmethod
    def write_to_scistree(genotype_matrix):
        """
        IO: Write the genotype matrix to a file.
        """
        nsite, ncell = genotype_matrix.shape
        prefix = uuid.uuid4()
        output = f"{prefix}.scistree.out"
        with open(output, "w") as out:
            out.write(f"HAPLOID\n")
            for i in range(nsite):
                for j in range(ncell):
                    prob = genotype_matrix[i, j]
                    out.write(f" {prob:.5f}")
                out.write("\n")
        return output

    @staticmethod
    def read_scistree_genotype(prefix):
        """
        Get the genotype matrix from the outputs of Scistree2.
        """
        geno_file = f"{prefix}.genos.imp"
        genotypes = []
        with open(geno_file, "r") as f:
            for line in f.readlines():
                if line.startswith("Site"):
                    line = line.strip()
                    genos = line.split("\t")[1].split()
                    genos = list(map(int, genos))
                    genotypes.append(genos)
        return np.array(genotypes)

    def bootstrap(self, tree, gp, num_bootstrap=100, num_site=-1):
        """
        Felsenstein’s tree bootstrp.
        """
        if num_site == -1:
            num_site = gp.nsite
        clades = {}
        for i in range(num_bootstrap):
            subgp = gp.subsample(n=num_site)
            t, geno, ml = self.infer(subgp)
            for split in t.get_splits():
                if split in clades:
                    clades[split] += 1
                else:
                    clades[split] = 1
        for node in tree.get_all_nodes():
            if tree[node].is_leaf() or tree[node].is_root():
                tree[node].branch_confidence = 1.0
            else:
                c = frozenset([leaf.name for leaf in tree[node].get_leaves()])
                if c in clades:
                    n = clades[c]
                else:
                    n = 0
                tree[node].branch_confidence = n / num_bootstrap
        return tree

    def infer(self, gp, verbose=False):
        """
        Run Scistree2 local search.

        Args:
            gp [scistree2.probability.GenotypeProbability]: Genotype probability matrix.
            verbose [bool]: Show logs.

        Returns:
            tree [scistree2.tree.BaseTree]: The optimal tree.
            imputed_genotype [pandas.DataFrame]: Imputed genotype. 0 as wild type, 1 as mutation.
            ml [float]: Log Likelihood of the optimal tree.
        """
        cell_names = gp.cell_names
        output = self.write_to_scistree(gp.probs)
        cmd = self.cmd + [f"{output}"]
        cmd = " ".join(cmd)
        try:
            res = (
                sp.run(cmd, shell=True, stdout=sp.PIPE, encoding="utf-8")
                .stdout.strip()
                .split("\n")
            )
            if verbose:
                print("\n".join(res))
            nwk = res[-2].split(":")[1].strip() + ";"
            if cell_names:
                t = relabel(
                    from_newick(nwk),
                    name_map={str(i + 1): name for i, name in enumerate(cell_names)},
                )
                nwk = t.output()
            imp_geno, ml, tree = evaluate(gp, nwk, return_tree=True)
            return tree, imp_geno, ml
        # clean output
        except Exception as e:
            print("scistree running failed.")
            if os.path.exists(output):
                os.remove(output)
            raise e
        finally:
            if os.path.exists(output):
                os.remove(output)
            if os.path.exists(f"{output}.genos.imp"):
                os.remove(f"{output}.genos.imp")


@staticmethod
def evaluate(gp, nwk, return_tree=False):
    """
    Evaluate a tree given genotype probabilities and return the imputated genotype.

    Args:
        gp [scistree2.probablity.GenotypeProbablity]: Genotype probability matrix.
        nwk [str]: Newick string of the optimal tree.
        return_tree: If return the optimal tree.

    Returns:
        imputed_genotype [pandas.DataFrame]: Imputed genotype. 0 as wild type, 1 as mutation.
        ml [float]: Log likelihood of the optimal tree.
        tree: Optinal
    """
    assert isinstance(nwk, str), "tree should be a newick string."
    cell_names = gp.cell_names
    site_names = gp.site_names
    tree = relabel(
        from_newick(nwk), name_map={name: str(i) for i, name in enumerate(cell_names)}
    )
    traveror = TraversalGenerator(order="post")
    max_mls = np.zeros(gp.probs.shape[0])  # - np.inf
    max_ml_nodes = [None] * gp.probs.shape[0]
    g = np.log(1 - gp.probs) - np.log(
        gp.probs
    )  # in log space to avoid numerical overflow
    for node in traveror(tree):
        if node.is_leaf():
            likelihood = g[:, int(node.name)]
        else:
            likelihood = (
                node.get_children()[0].likelihood + node.get_children()[1].likelihood
            )
        for i, l in enumerate(likelihood):
            if l > max_mls[i]:
                max_mls[i] = l
                max_ml_nodes[i] = node
        node.likelihood = likelihood
    max_mls += np.log(gp.probs).sum(axis=1)
    imputed_genotype = np.zeros_like(gp.probs, dtype=int)
    # print(imputed_genotype.sum(axis=-1))
    for i, ml_node in enumerate(max_ml_nodes):
        if ml_node:
            ml_node.add_mutation(site_names[i])
            inds = [int(leaf.name) for leaf in ml_node.get_leaves()]
            imputed_genotype[i, inds] = 1
    # rename leaves back to cell names
    tree = relabel(tree, name_map={str(i): name for i, name in enumerate(cell_names)})
    if return_tree:
        return (
            pd.DataFrame(imputed_genotype, index=gp.site_names, columns=gp.cell_names),
            sum(max_mls[max_mls != -np.inf]),
            tree,
        )
    return pd.DataFrame(
        imputed_genotype, index=gp.site_names, columns=gp.cell_names
    ), sum(max_mls[max_mls != -np.inf])
