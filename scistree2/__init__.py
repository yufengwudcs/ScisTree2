import os
# from ... import popgen
import numpy as np 
from . import metric
from . import probability
from . import treeutils as util
from .scistree import ScisTree2, evaluate
from .reader import read_vcf



### version updates ###
# __version__ = "0.1.0" # original submission
# __version__ = "0.2.0" # add ete3 tree, add mutation profile at branches.
# __version__ = "0.3.0" # add bootstrapping.
# __version__ = "0.4.0" # to pypi.
__version__ = "0.5.0" # add win support.