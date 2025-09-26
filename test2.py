from phytreeviz import TreeViz, load_example_tree_file
import popgen

nwk = '(((1)2,2)2:3,(3,4)4:2)5:1;'
# nwk = popgen.utils.get_random_binary_tree(100).output()
tv = TreeViz(nwk)
tv.show_branch_length()
tv.show_confidence()
tv.savefig("api_example01.png", dpi=300)