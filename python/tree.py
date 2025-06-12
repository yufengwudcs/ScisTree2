from collections import OrderedDict, defaultdict
from .node import Node
import pptree
import copy
import _pickle as cPickle


def deepcopy(obj):
    try:
        return cPickle.loads(cPickle.dumps(obj, -1))
    except Exception:
        return copy.deepcopy(obj)

class BaseTree(object):
    """
    A tree class, for extracting features from genealogical trees in the future.
    Arguments:
        args:
    # >>> tree = BaseTree();
    # >>> tree = popgen.utils.treeutils.from_newick('((1,2),(3,4));') # or
    """
    def __init__(self, root=None):
        self.root = root
        self._update()


    def __contains__(self, node):
        # if a node is in this tree
        if isinstance(node, Node):
            node = node.identifier
        if node in self._nodes:
            return True
        else:
            return False

    def __getitem__(self, identifier) -> Node:
        # get a node from node list
        return self._nodes[identifier]

    def __len__(self):
        # obtain number of nodes in a tree
        return len(self._nodes)
    
    def __eq__(self, tree):
        anchor = tree[tree.get_leaves()[0]].name
        return tree.get_splits(True, anchor) == self.get_splits(True, anchor)

    def __hash__(self):
        return hash(self.get_splits(True))

    def __str__(self):
        return self.output()

    def _update(self):
        if self.root is not None:
            self._nodes = self.root.get_descendants_dict()
            self._nodes[self.root.identifier] = self.root


    def create_node(self, identifier=None, name=None, parent=None) -> Node:
        node = Node(identifier=identifier, name=name)
        self.add_node(node, parent)
        return node

    def add_node(self, node, parent=None):
        if node.identifier in self._nodes:
            raise Exception('cannot add the node that has already been in the tree.')
        if parent is None:
            if self.root:
                raise Exception('root has already existed and parent cannot be none.')
            self.root = node
            self._nodes[node.identifier] = node
            # /* set root and make its level as 0 */
            node.set_parent(None)
            return
        pid = parent.identifier if isinstance(parent, Node) else parent
        if not pid in self._nodes:
            raise Exception('parent not found in this tree.')
        self._nodes[node.identifier] = node
        # /* link node with its parent */
        node.set_parent(self[pid])
        self[pid].add_child(node)
        return

    def get_all_nodes(self):
        return self._nodes

    def get_leaves(self, return_label=False):
        leaves = [node.name if return_label else node.identifier for node in self.root.get_leaves()]
        return leaves

    def get_splits(self, return_label=False, contains_leaf=None):
        def get_complement_split(split):
            return frozenset(filter(lambda x: x not in split, leaves))
        leaves = self.get_leaves(return_label=True)
        splits = set()
        for nid in self.get_all_nodes():
            # if it is root or a leaf which means trivial split, pass
            if not self._nodes[nid].is_leaf() and not self._nodes[nid].is_root():
                split = frozenset([node.name if return_label else node.identifier for node in self._nodes[nid].get_leaves()])
                if 1 < len(split) < len(leaves)-1 :
                    if contains_leaf and contains_leaf not in split:
                        splits.add(get_complement_split(split))
                    else:
                        splits.add(split)
        return frozenset(splits)

    def to_dict(self):
        # return a dict for the whole tree.
        pass

  
    def output(self, output_format='newick_sorted', branch_lengths=False):
        def _newick_unsorted(node, branch_lengths):
            if node.is_leaf():
                return node.name
            fstr = '(' + ','.join(['{newick}:{branch}'
                                  .format(newick=_newick_unsorted(child, branch_lengths), branch=child.branch)
                                   if branch_lengths else _newick_unsorted(child, branch_lengths)
                                   for child in node.get_children()]) + ')'
            return fstr

        def _newick_sorted(node, branch_lengths):
            if node.is_leaf():
                return node.name
            newick_children = []
            for child in node.get_children():
                if branch_lengths:
                    newick_children.append('{newick}:{branch}'.format(newick=_newick_sorted(child, branch_lengths), branch=child.branch))
                else:
                    newick_children.append(_newick_sorted(child, branch_lengths))
            fstr = '(' + ','.join(sorted(newick_children)) + ')'
            return fstr
        
        def newick():
            return _newick_unsorted(self.root, branch_lengths) + ';'
        
        def newick_sorted():
            return _newick_sorted(self.root, branch_lengths) + ';'

        funcs = {'newick': newick, 'newick_sorted': newick_sorted}
        return funcs[output_format]()

    def draw(self, **kwargs):
        pptree.print_tree(self.root, "_children", **kwargs)

    def copy(self):
        return deepcopy(self)