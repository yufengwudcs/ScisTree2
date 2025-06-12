# Here, we provide a Python interface for conveniently calling Scistree2.
import os
import multiprocessing as mp
import uuid
import numpy as np
import subprocess as sp
from .treeutils import *
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

class ScisTree2():
    def __init__(self,
                 threads=-1, nj=False, spr=True, nni=False, iterative=False, verbose=True):
        self.bin_path = os.path.join(os.path.dirname(__file__), 'bin', 'scistree')
        self.nj = nj
        self.spr = spr
        self.nni = nni
        self.cmd = self.build_cmd(self.bin_path, threads, nj, nni, verbose)
        
    def build_cmd(self, bin_path, threads, nj, nni, verbose):
        if threads == -1:
            threads = mp.cpu_count()
        cmd = [bin_path, '-T', str(threads)]
        if verbose:
            cmd.append('-v')
        if nj:
            cmd.append('-n')
        if nni:
            cmd.append('-q')
        return cmd 
    
    """
    IO: Write the genotype matrix to a file.
    """
    @staticmethod
    def write_to_scistree(genotype_matrix):
        nsite, ncell = genotype_matrix.shape
        prefix = uuid.uuid4()
        output = f'{prefix}.scistree.out'
        with open(output, 'w') as out:
            out.write(f'HAPLOID {nsite} {ncell}')
            for i in range(ncell):
                out.write(f' {i}')
            out.write('\n')
            for i in range(nsite):
                out.write(f's{i}')
                for j in range(ncell):
                    prob = genotype_matrix[i, j]
                    out.write(f' {prob:.10f}')
                out.write('\n')
        return output
    
    """
    Get the genotype matrix from the outputs of Scistree2.
    """
    @staticmethod
    def read_scistree_genotype(prefix):
        geno_file = f'{prefix}.genos.imp'
        genotypes = []
        with open(geno_file, 'r') as f:
            for line in f.readlines():
                if line.startswith('Site'):
                    line = line.strip()
                    genos = line.split('\t')[1].split()
                    genos = list(map(int, genos))
                    genotypes.append(genos)
        return np.array(genotypes)
    
    """
    Run Scistree2 local search.

    Args:
        geno [numpy.ndarray]: Genotype probability matrix.

    Returns:
        imputed_genotype [numpy.ndarry]: Imputed genotype. 0 as wild type, 1 as mutation.
        nwk [str]: Newick string of the optimal tree.
        ml [float]: Log Likelihood of the optimal tree.
        time [int]: Running time in seconds.
    """
    def infer(self, geno):
        output = self.write_to_scistree(geno)
        cmd = self.cmd + [f'{output}']
        cmd = ' '.join(cmd)
        try:
            res = sp.run(cmd, shell=True, stdout=sp.PIPE, encoding='utf-8').stdout.strip().split('\n')
            if not self.nj:
                # time = float(res[-1].split('=')[1].split('seconds')[0])
                nwk = res[-2].split(':')[1].strip() + ';'
                ml = float(res[-4].split(',')[0].split(':')[1])
                imp_geno = self.read_scistree_genotype(output)
                return imp_geno, nwk, ml
            else:
                # time = float(res[-1].split('=')[1].split('seconds')[0])
                nwk = res[-2].split(':')[1].strip() + ';'
                return nwk
        except Exception as e:
                print('scistree running failed.')
                raise e
        finally:
            if os.path.exists(output):
                os.remove(output)
            if os.path.exists(f'{output}.genos.imp'):
                os.remove(f'{output}.genos.imp')
    """
    Evaluate a tree given genotype probabilities and return the imputated genotype.

    Args:
        geno [numpy.ndarray]: Genotype probability matrix.
        nwk [str]: Newick string of the optimal tree.

    Returns:
        imputed_genotype [numpy.ndarry]: Imputed genotype. 0 as wild type, 1 as mutation.
        ml [float]: Log likelihood of the optimal tree.
    """
    @staticmethod
    def evaluate(geno, nwk, offset=-1):
        assert isinstance(nwk, str), 'tree should be a newick string.'
        tree = relabel(from_newick(nwk), offset=offset)
        traveror = TraversalGenerator(order='post')
        max_mls = np.zeros(geno.shape[0]) # - np.inf 
        max_ml_nodes = [None] * geno.shape[0]
        g = np.log(1-geno) - np.log(geno) # in log space to avoid numerical overflow
        for node in traveror(tree):
            if node.is_leaf():
                likelihood = g[:, int(node.name)]
            else:
                likelihood = node.get_children()[0].likelihood + node.get_children()[1].likelihood
            for i, l in enumerate(likelihood):
                if l > max_mls[i]:
                    max_mls[i] = l
                    max_ml_nodes[i] = node
            node.likelihood = likelihood
            # print(likelihood)
        max_mls += np.log(geno).sum(axis=1)
        imputed_genotype = np.zeros_like(geno, dtype=int)
        # print(imputed_genotype.sum(axis=-1))
        for i, ml_node in enumerate(max_ml_nodes):
            if ml_node:
                inds = [int(leaf.name) for leaf in ml_node.get_leaves()]
                imputed_genotype[i, inds] = 1
        return imputed_genotype, sum(max_mls[max_mls != -np.inf])

