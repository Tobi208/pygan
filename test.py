import unittest
# from algorithms.lca_original import *
from algorithms.lca import *
from graph.graph import *
from graph.node import Node
from graph.edge import Edge


def generate_graph(height):
    num_nodes = 2 ** (height + 1) - 1
    g = Graph()
    g.nodes = [Node(i) for i in range(num_nodes)]
    edges = [[Edge(i * 2, g.nodes[i], g.nodes[i * 2 + 1]), Edge(i * 2 + 1, g.nodes[i], g.nodes[i * 2 + 2])]
             for i in range(2 ** height - 1)]
    g.edges = [item for sublist in edges for item in sublist]
    for e in g.edges:
        e.source.out_edges.append(e)
        e.target.in_edges.append(e)
    g.root = g.nodes[0]
    return g


g = generate_graph(4)
n = g.nodes


class PyganTests(unittest.TestCase):

    def test_compute_addresses(self):

        id2a = {}
        a2id = {}
        compute_addresses(g, id2a, a2id)

        for k, v in id2a.items():
            self.assertEqual(a2id[v], k)
        #
        # self.assertEqual(list(map(ord, id2a[0])), [])
        # self.assertEqual(list(map(ord, id2a[1])), [1])
        # self.assertEqual(list(map(ord, id2a[14])), [2, 2, 2])
        # self.assertEqual(list(map(ord, id2a[15])), [1, 1, 1, 1])
        # self.assertEqual(list(map(ord, id2a[17])), [1, 1, 2, 1])

        self.assertEqual(id2a[0], ())
        self.assertEqual(id2a[1], (0,))
        self.assertEqual(id2a[14], (1, 1, 1))
        self.assertEqual(id2a[15], (0, 0, 0, 0))
        self.assertEqual(id2a[17], (0, 0, 1, 0))

    def test_lcs_addressing(self):

        id2a = {}
        a2id = {}
        compute_addresses(g, id2a, a2id)

        # lca of root + other node
        # returns other node instead of root in original algorithm
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[0], id2a[15]], ignore_ancestors=True)]).node_id, 15)
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[0], id2a[15]])]).node_id, 0)
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[0], id2a[15], id2a[18]])]).node_id, 0)

        self.assertIs(g.get_node(a2id[get_common_prefix([])]).node_id, 0)
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[0]])]).node_id, 0)
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[0], id2a[0]])]).node_id, 0)
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[16]])]).node_id, 16)
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[1], id2a[3]])]).node_id, 1)
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[11], id2a[12]])]).node_id, 5)
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[17], id2a[19]])]).node_id, 1)
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[14], id2a[15]])]).node_id, 0)
        self.assertIs(g.get_node(a2id[get_common_prefix([id2a[15], id2a[16], id2a[17]])]).node_id, 3)


if __name__ == '__main__':
    unittest.main()
