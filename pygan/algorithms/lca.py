from typing import List, Tuple, Dict

from pygan.datastructure.directed_graph import DirectedGraph, DirectedNode


def compute_addresses(tree: DirectedGraph, id2address: Dict[int, Tuple], address2id: Dict[Tuple, int]):
    """
    Computes node addresses used to compute LCA

    :param tree: phylotree
    :param id2address: map of ids to addresses
    :param address2id: map of addresses to ids
    """
    root = tree.root
    if root:
        build_id2address_rec(root, (), id2address, address2id)


def build_id2address_rec(v: DirectedNode, path: Tuple, id2address: Dict[int, Tuple], address2id: Dict[Tuple, int]):
    """
    Computes the id to address mapping and vice versa

    :param v: current node
    :param path: sequence of junctions taken to reach node
    :param id2address: map of ids to addresses
    :param address2id: map of addresses to ids
    """
    id2address[v.node_id] = path
    address2id[path] = v.node_id
    for i, e in enumerate(v.out_edges):
        build_id2address_rec(e.target, path + (i,), id2address, address2id)


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

# def get_common_prefix_alt(addresses, ignore_ancestors=False):
#
#     import itertools
#
#     addresses = [address for address in addresses if address is not None]
#
#     if len(addresses) == 0:
#         return ()
#     elif len(addresses) == 1:
#         return addresses[0]
#
#     if ignore_ancestors:
#         addresses_zipped = itertools.zip_longest(*addresses)
#         for i, ps in enumerate(addresses_zipped):
#             non_none = [j for j in range(len(ps)) if ps[j] is not None]
#             if len(non_none) == 1:
#                 return addresses[non_none[0]]
#             if any(map(lambda p: p != ps[non_none[0]], [ps[j] for j in non_none])):
#                 return addresses[non_none[0]][:1]
#     else:
#         enum = enumerate(zip(*addresses))
#         for i, ps in enum:
#             if any(map(lambda p: p != ps[0], ps)):
#                 return addresses[0][:i]
#         return addresses[0][:len(list(enum)) + 1]
#     return ()
