from .node import Node 
from .tree import BaseTree
import warnings
import time
import random
import numpy as np
import os 
import subprocess as sp
"""
author: Haotian
created_at: 12/4/2021
description: utils for trees
"""

class TraversalGenerator(object):
    """
        A generator class used for tree traversal
        Arguments:
            order: traversal order
        >>> generator = TraversalGenerator(order='post')
        >>> for node in generator(tree):
                # to do something
                pass
    """

    def __init__(self, order='post'):
        self.order = order
        self.iterator = None

    def __call__(self, tree, order='post'):
        # calling function
        valid_methods = {'pre': self._pre, 'in': self._in, 'post': self._post}
        assert order in valid_methods, "order should be in ['pre', 'mid', 'post']"
        method = valid_methods[order]
        iterator = iter(method(tree))
        return iterator

    def _pre(self, tree):
        # pre-order traverse
        node = tree.root
        traverse_nodes = []
        while len(traverse_nodes) != len(tree):
            while not node.is_leaf():
                if node.identifier not in traverse_nodes:
                    traverse_nodes.append(node.identifier)
                    yield node
                flag = False
                for child in node.get_children():
                    if child.identifier not in traverse_nodes:
                        flag = True
                        node = child
                        break
                if not flag: node = node.parent
            if node.identifier not in traverse_nodes:
                traverse_nodes.append(node.identifier)
                yield node
            node = node.parent
                

    def _in(self, tree):
        # mid-order traverse
        warnings.warn("no implementation currently")

    def _post(self, tree):
        node = tree.root
        traverse_nodes = []
        while len(traverse_nodes) != len(tree):
            while node.get_children():
                flag = True
                for child in node.get_children():
                    if child.identifier not in traverse_nodes:
                        node = child
                        flag = False
                        break
                if flag:
                    break
            traverse_nodes.append(node.identifier)
            yield node
            node = node.parent

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self, method):
        valid_methods = {'pre': self._pre, 'in': self._in, 'post': self._post}
        if method in valid_methods:
            self._order = method
            self._method = valid_methods[method]
        else:
            raise Exception("order should be in ['pre', 'mid', 'post']")


def from_node(node: Node) -> BaseTree:
    """ Build a tree by directly setting a root """
    tree = BaseTree(root=node)
    return tree


def from_newick(newick: str) -> BaseTree:
    """ Build a tree according to a newick-format string """
    def _isvalid(s):
        checking_stack = []
        for ch in s:
            checking_stack.append(ch) if ch == '(' else None
            if ch == ')':
                if checking_stack:
                    checking_stack.pop()
                else:
                    return False
        return True if not checking_stack and ch == ';' else False

    def _next(i):
        stop_words = [',', ')', ';']
        if newick[i] in stop_words:
            return None, 0, i
        br = 0
        j = i + 1
        while newick[j] not in stop_words:
            j += 1
        if ':' in newick[i:j]:
            nid, br = newick[i:j].split(':')
        else:
            nid = newick[i:j]
        return nid.strip(), br, j

    newick = newick.strip()
    assert isinstance(newick, str), Exception('newick should be a string.')
    assert _isvalid(newick), Exception('invalid newick string.')
    nodes = []
    level, i, key = 0, 0, ''
    while not newick[i] == ';':
        if newick[i] in [',', ' ']:
            i += 1
            continue
        if newick[i] == '(':
            level += 1
            i += 1
            continue
        if newick[i] == ')':
            if not nodes:
                raise Exception('newick: bad parsing.')
            identifier, branch, end = _next(i + 1)
            node = Node(identifier=identifier, branch=branch)
            while nodes:
                if nodes[-1][1] == level:
                    child = nodes.pop()
                    node.add_child(child[0])
                    child[0].set_parent(node)
                else:
                    break
            level -= 1
            nodes.append((node, level))
            i = end
            continue
        identifier, branch, end = _next(i)
        node = Node(identifier=identifier, branch=branch)
        nodes.append((node, level))
        i = end

    root = nodes.pop()
    tree = from_node(root[0])
    return tree


# perturb the branch lengths for data without clock property
def perturb_tree_length(tree, min=.5, max=1.5, mode='mul'):
    tree = tree.copy()
    for node_str in tree.get_all_nodes():
        node = tree[node_str]
        if not node.is_root():
            try:
                if mode == 'mul':
                    rand = np.random.uniform(min, max)
                    node.branch = np.round(node.branch * rand, 4)
                if mode == 'add':
                    rand = np.random.rand()
                    node.branch = np.round(node.branch + rand, 4)
            except ValueError:
                raise Exception('branch is not a valid number')
    return tree


# relabel the leaves, this is usually used when there are different start indices
# only be using when labels are integers
# name_map should be a dictionary which is formatted as {'old_name': 'new_name'}
def relabel(tree, offset=0, name_map=None):
    tree = tree.copy()
    # try:
    if offset != 0 and name_map:
        raise Exception('error, both offset and name_map are specified')
    elif offset != 0:
        for leaf in tree.get_leaves():
            tree[leaf].name = str(int(tree[leaf].name)+offset)
    elif offset == 0 and name_map == None:
        pass
    else:
        for leaf in tree.get_leaves():
            name = tree[leaf].name
            if name in name_map:
                tree[leaf].name = name_map[name]
    # except Exception as e:
    #     print(e)
    #     exit()
    return tree


"""
    Construct a perfect phylogeny (no conflicts) (Dan Gusfield/ https://doi.org/10.1002/net.3230210104)
    Author: haotian Z
    Date: 12/27/22
"""
class BNode(Node):
    # extended Node class for recording mutations 
    def __init__(self, identifier=None, name=None, branch=0):
        super().__init__(identifier, name, branch)
        self.mutations = []

    def add_mutations(self, mutations):
        self.mutations.extend(mutations)


#  vanilla implementation takes O(nm^2) that simply check every pair of columns to see if they are disjoint or one includes the other one.
def check_no_conflict_vanilla(mat):
    nrow, ncol = mat.shape
    gametes = [[0,1], [1,0], [1,1]]
    for i in range(ncol):
        for j in range(ncol):
            flag = False
            for gamete in gametes:
                # print(gamete, mat[:, [i,j]])
                if gamete not in mat[:, [i,j]].tolist():
                    flag = True
                    break
            if not flag:
                return False
    return True


def remove_homozygous_columns(mat):
    mask = mat.sum(axis=0) > 0
    return mat[:, mask]


# useless one
def binary_number(arr):
    s = 0
    for i, ele in enumerate(arr[::-1]):
        if ele:
            s += 2**i
    return s 


# sort by binary numbers and remove duplicates
def rearrangement(mat):
    nrow, ncol = mat.shape
    binary_strs = []
    for i in range(ncol):
        binary_strs.append(''.join([str(val) for val in mat[:, i]]))
    sorted_index = sorted(range(ncol), key=lambda x:binary_strs[x], reverse=True)
    # build index and remove duplicates
    groups = {}
    group = []
    final_index = []
    for i in range(ncol):
        group.append(sorted_index[i])
        if i==ncol-1 or binary_strs[sorted_index[i]] != binary_strs[sorted_index[i+1]]:
            final_index.append(sorted_index[i])
            groups[sorted_index[i]] = group
            group = []
    return final_index, groups
    

# get a matrix for determining L(j)
def preprocess(mat):
    nrow, ncol = mat.shape
    indices, groups = rearrangement(mat)
    mat_ = mat[:, indices].copy()
    # for loop to check every entity with value 1
    pre_mat = np.zeros_like(mat_)
    for i in range(mat_.shape[0]):
        pre = 0
        for j in range(mat_.shape[1]):
            if mat_[i,j] == 1:
                pre_mat[i,j] = pre
                pre = j+1
    return mat_, pre_mat, indices, groups


#  Gusfield here did some sortings based on binary numbers then make an improvement to O(nm) [n: # of cells, m: # of columns]
def check_no_conflict(mat):
    nrow, ncol = mat.shape
    mat_, L, indices, groups = preprocess(mat)
    for j in range(L.shape[1]):
        maxx = L[:, j].max()
        for i in range(L.shape[0]):
            if mat_[i,j] == 1 and L[i,j] != maxx:
                return False
    return True
                    

# phylogeny construction from Gusfield's paper from a conflict-free matrix
def build_perfect_phylogeny(mat):
    mat = remove_homozygous_columns(mat)
    if not check_no_conflict(mat):
        raise Exception('This are some conflicts in the matrix provided.')
    mat_, L, indices, groups = preprocess(mat)
    L = np.max(L, axis=0)
    # step 1: build mutation tree
    root = BNode(identifier='root')
    mutation_nodes = {}
    for mutation in groups:
        node = BNode()
        # node = BNode(name=str(groups[mutation]))
        node.add_mutations(groups[mutation])
        mutation_nodes[mutation] = node
    for j in range(len(L)):
        if L[j] == 0:
            mutation = indices[j]
            root.add_child(mutation_nodes[mutation])
            mutation_nodes[mutation].set_parent(root)
        else:
            mutation_from = indices[L[j]-1]
            mutation_to = indices[j]
            node_from = mutation_nodes[mutation_from]
            node_to = mutation_nodes[mutation_to]
            node_from.add_child(node_to)
            node_to.set_parent(node_from)

    # step 2: add cells into the mutation tree
    nrow, ncol = mat_.shape
    for i in range(nrow):
        max_index = ncol-np.argmax(mat_[i][::-1])-1
        last_mutation_node = mutation_nodes[indices[max_index]]
        leaf_node = BNode(identifier=i)
        last_mutation_node.add_child(leaf_node)
        leaf_node.set_parent(last_mutation_node)
    
    # phylogeny = popgen.utils.from_node(root)
    # phylogeny.print()
    # step 3: tree prune
    for label in mutation_nodes:
        node = mutation_nodes[label]
        if len(node.get_children()) == 1:
            child = node.get_children()[0]
            node.parent.children.pop(node.identifier)
            node.parent._children.remove(node)
            node.parent.add_child(child)
            child.set_parent(node.parent)
            child.add_mutations(node.mutations)
    phylogeny = from_node(root)
    return phylogeny

def single_spr_move(tree, verbose=False):
    """
        Note: only for binary tree
        a spr move is to randomly select a node and swap its parent with another node
        the parent node should not be the root and the other node cannot be its descendant
        return a new tree

        1. allow the subtree directly under the root
    """
    
    while True:
        tree = tree.copy()
        nodes = list(tree.get_all_nodes()) # identifiers
        node = tree[random.choice(nodes)]
        parent = node.parent
        if node.is_root():
            continue
        sibling = node.get_siblings()[0]
        if not parent.is_root():
            grandparent = parent.parent
            grandparent.remove_child(parent)
            parent.set_parent(None)
            parent.remove_child(sibling)
            sibling.set_parent(None)
            grandparent.add_child(sibling)
            sibling.set_parent(grandparent)
        else:
            parent.remove_child(sibling)
            sibling.set_parent(None)
            tree = BaseTree(sibling)
        tree._update()
        candidates = list(tree.get_all_nodes().values())
        if len(candidates) == 1:
            continue
        break
    if verbose:
        print('spr remove: ', node.identifier)
    while True:
        new_sibling = random.choice(candidates)
        if new_sibling == sibling:
            continue
        if verbose:
            print('spr regraft: ', new_sibling.identifier)
        if new_sibling.is_root():
            parent.add_child(new_sibling)
            new_sibling.set_parent(parent)
            tree = BaseTree(parent)
        else:
            parent.add_child(new_sibling)
            new_sibling.parent.remove_child(new_sibling)
            new_sibling.parent.add_child(parent)
            parent.set_parent(new_sibling.parent)
            new_sibling.set_parent(parent)        
        break
    tree._update()
    return tree

# def single_spr_move(tree):
#     """
#         Note: only for binary tree
#         a spr move is to randomly select a node and swap its parent with another node
#         the parent node should not be the root and the other node cannot be its descendant
#         return a new tree
#     """
#     tree = tree.copy()
#     nodes = list(tree.get_all_nodes()) # identifiers
    
#     while True:
#         node = tree[random.choice(nodes)]
#         parent = node.parent
#         if node.is_root():
#             continue
#         if parent.is_root():
#             continue
#         # get all nodes that are not descendants of node
#         candidates = []
#         for idx in nodes:
#             n = tree[idx]
#             if n not in node.get_descendants() and n != node:
#                 candidates.append(n)
#         if not candidates:
#             continue
                
#         new_sibling = random.choice(candidates)
#         if new_sibling.is_root():
#             continue
#         if new_sibling == parent:
#             continue
#         # swap
#         # print('node', node.identifier)
#         # print('sib', new_sibling.identifier)  
#         grandparent = parent.parent
#         sibling = parent.get_children()[0]
#         grandparent.remove_child(parent)
#         parent.set_parent(None)
#         parent.remove_child(sibling)
#         grandparent.add_child(sibling)
#         sibling.set_parent(grandparent)
#         sibling_parent = new_sibling.parent
#         sibling_parent.remove_child(new_sibling)
#         new_sibling.set_parent(None)
#         parent.add_child(new_sibling)
#         new_sibling.set_parent(parent)
#         sibling_parent.add_child(parent)
#         parent.set_parent(sibling_parent)
#         break

#     tree._update()
#     return tree


def spr_move(tree, move):
    for _ in range(move):
        tree = single_spr_move(tree)
    return tree


def spr_distance_cpp(newick1, newick2, upper_bound=None, spr_exec_path='/home/haz19024/softwares/rspr/rspr'):
    """
    Computes the SPR distance between two vectors.
    """
    assert os.path.exists(spr_exec_path), "rSPR binary file not found."
    if not upper_bound:
        a = sp.run([spr_exec_path, "-bb", "-q"], input=f'{newick1}\n{newick2}'.encode(), stdout=sp.PIPE)
    else:
        a = sp.run([spr_exec_path, "-bb", "-q", "-split-approx", str(upper_bound)], input=f'{newick1}\n{newick2}'.encode(), stdout=sp.PIPE)
    result = a.stdout.decode().split('\n')[-2].split('=')[1]
    return float(result) 


def spr_distance(tree1, tree2):    
    newick1 = str(tree1)
    newick2 = str(tree2)
    return spr_distance_cpp(newick1, newick2)

def get_random_binary_tree(n_leave, start_index=0, random_branch=True):
    # create a random binary tree given the number of leaves
    nodes = [BNode(identifier=i) for i in range(start_index, n_leave+start_index)]
    while len(nodes) > 1:
        i = np.random.randint(0, len(nodes))
        j = np.random.randint(0, len(nodes))
        while j == i:
            j = np.random.randint(0, len(nodes))
        node = Node()
        node.add_child(nodes[i])
        node.add_child(nodes[j])
        if random_branch:
            nodes[i].branch = np.round(np.random.rand(), 4)
            nodes[j].branch = np.round(np.random.rand(), 4)
        nodes[i].set_parent(node)
        nodes[j].set_parent(node)
        nodes.pop(max(i,j))
        nodes.pop(min(i,j))
        nodes.append(node)
    root = nodes.pop()
    tree = from_node(root)
    return tree


def build_no_repeat_clade(g):
    clades = set()


def mutated_clade_dist(g1, g2):
    pass