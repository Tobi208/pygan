from pygan.tree.phylo_tree import PhyloTree, PhyloNode


def get_phylo_tree(file: str) -> PhyloTree:
    """
    Read a file containtng a newick tree in format (1,2,(3,4)5)6  and parse it to phylogenetic tree

    :param file: file containing newick tree in format (1,2,(3,4)5)6
    :return: phylogenetic tree
    """
    return parse(read(file))


def read(file: str) -> str:
    """
    Read a file

    :param file: filepath
    :return: content of the file
    """
    with open(file, 'r') as f:
        return f.read()


def parse(newick: str) -> PhyloTree:
    """
    Parse a newick string to a phylogentic tree iteratively in O(n)
    The newick tree is trusted to be in the format (1,2,(3,4)5)6

    May be extended to include depths of nodes

    :param newick: newick tree in format (1,2,(3,4)5)6
    :return: phylogenetic tree
    """

    # adjust parent/node pointers per bracket
    tree = PhyloTree()
    parent = None
    node = None
    tax_id = ''

    for c in newick:

        # most likely case to happen, check first
        # just append the id digit to the tax_id
        if c.isdigit():

            tax_id += c

        # , indicates a sibling, don't change level
        # node is left of ,
        elif c == ',':

            # check if node was already created in a preceding iteration
            # if not, create it
            if not node:
                node = PhyloNode()
                node.parent = parent

            # add entry to parent and tree
            tax_id = int(tax_id)
            node.tax_id = tax_id
            tree.nodes[tax_id] = node
            parent.children.append(node)

            # flush id and node
            tax_id = ''
            node = None

        # ) indicates conclusion of a child, go back to a higher level
        # node is left of )
        elif c == ')':

            # check if node was already created in a preceding iteration
            # if not, create it
            if not node:
                node = PhyloNode()
                node.parent = parent

            # add entry to parent and tree
            tax_id = int(tax_id)
            node.tax_id = tax_id
            tree.nodes[tax_id] = node
            parent.children.append(node)

            # flush id and update node/parent pointers to go one level up
            tax_id = ''
            node = parent
            parent = node.parent

        # ( indicates creation of a child, go to a lower level
        elif c == '(':

            # node is definitely new
            node = PhyloNode()

            # if parent exists (root.parent = None), set pointer
            if parent:
                node.parent = parent

            # update node/parent pointers to go one level down
            parent = node
            node = None

    # first node created, but last node concluded is root
    tax_id = int(tax_id)
    node.tax_id = tax_id
    tree.nodes[tax_id] = node
    tree.root = node

    return tree


def to_newick(node: PhyloNode) -> str:
    """
    Parse a phylogenetic tree back to a newick tree in format (1,2,(3,4)5)6.
    Call with root node.

    Mainly for testing/debugging. Implementation is highly inefficient.

    :param node: phylogentic tree node
    :return: newick tree in format (1,2,(3,4)5)6
    """
    if len(node.children) == 0:
        return str(node.tax_id)
    else:
        return '(' + ','.join([to_newick(child) for child in node.children]) + ')' + str(node.tax_id)