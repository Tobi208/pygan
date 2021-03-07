import sys


def compute_addresses(tree, id2address, address2id):
    root = tree.root
    if root:
        build_id2address_rec(root, '', id2address, address2id)


def build_id2address_rec(v, path, id2address, address2id):
    node_id = v.node_id
    id2address[node_id] = path
    address2id[path] = node_id
    if len(v.out_edges) < sys.maxunicode:
        count = 1
        for e in v.out_edges:
            build_id2address_rec(e.target, path + chr(count), id2address, address2id)
            count += 1
    else:
        count1 = 1
        count2 = 2
        for e in v.out_edges:
            if count1 == sys.maxunicode:
                count2 += 1
                count1 = 1
            build_id2address_rec(e.target, path + chr(count1) + chr(count2), id2address, address2id)
            count1 += 1


def get_common_prefix(addresses, ignore_ancestors):
    if len(addresses) == 0:
        return ''
    elif len(addresses) == 1:
        return addresses[0]

    reference = None

    for other in addresses:
        if other and len(other) > 0:
            if not reference:
                reference = other
            else:
                if ignore_ancestors and len(other) > len(reference) \
                        or not ignore_ancestors and len(other) < len(reference):
                    reference = other

    if not reference:
        return ''

    for (i, char) in enumerate(reference):
        for other in addresses:
            if other and i < len(other) and other[i] != char:
                return reference[:i]

    return reference
