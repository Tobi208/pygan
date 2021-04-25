from typing import List, Tuple, Dict

from pygan.tree.phylo_tree import PhyloTree, PhyloNode


def compute_addresses(tree: PhyloTree, id2address: Dict[int, Tuple], address2id: Dict[Tuple, int]):
    """
    Computes node addresses used to compute LCA

    :param tree: phylotree
    :param id2address: map of ids to addresses
    :param address2id: map of addresses to ids
    """
    root = tree.root
    if root:
        build_id2address_rec(root, (), id2address, address2id)


def build_id2address_rec(v: PhyloNode, path: Tuple, id2address: Dict[int, Tuple], address2id: Dict[Tuple, int]):
    """
    Computes the id to address mapping and vice versa

    :param v: current node
    :param path: sequence of junctions taken to reach node
    :param id2address: map of ids to addresses
    :param address2id: map of addresses to ids
    """
    id2address[v.tax_id] = path
    address2id[path] = v.tax_id
    for i, c in enumerate(v.children):
        build_id2address_rec(c, path + (i,), id2address, address2id)


def get_common_prefix(addresses: List[Tuple], ignore_ancestors: bool = False) -> Tuple:
    """
    Compute common prefix of a list of addresses

    :param addresses: list of addresses
    :param ignore_ancestors: use longest address or shortest address as reference
    :return: prefix
    """

    # filter out None types
    addresses = [address for address in addresses if address is not None]

    # catch length = 0 or 1
    if len(addresses) == 0:
        return ()
    elif len(addresses) == 1:
        return addresses[0]

    # if ancestors are to be considered, and root address is contained
    # return root address
    if not ignore_ancestors and any(map(lambda x: len(x) == 0, addresses)):
        return ()

    # if ancestors are to be ignored:
    #   use longest address as reference
    reference = addresses[0]
    if ignore_ancestors:
        for address in addresses:
            if len(address) > len(reference):
                reference = address
    # if ancestors are to be considered:
    #   use shortest address as reference
    else:
        for address in addresses:
            if len(address) < len(reference):
                reference = address

    # iterate until first difference in paths is discovered
    for (i, branch) in enumerate(reference):
        for address in addresses:
            if i < len(address) and address[i] != branch:
                return reference[:i]

    # if ancestors are to be ignored and all addresses are ancestors of reference
    # return reference
    return reference
