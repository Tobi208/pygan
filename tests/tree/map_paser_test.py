import unittest
import os
from pygan.tree.map_parser import *
from pygan.tree.newick_parser import parse as nwk_parse, get_phylo_tree


class MapParserTest(unittest.TestCase):

    def test_readlines(self):

        with self.assertRaises(FileNotFoundError):
            readlines('not a file.map')

        with open('test.map', 'w') as f:
            f.write('1\ta\n2\tb\n3\tc\n')
        self.assertListEqual(readlines('test.map'), [['1', 'a'], ['2', 'b'], ['3', 'c']])
        os.remove('test.map')

    def test_parse(self):

        tree = nwk_parse('(1,2,(3,4)5)6')
        lines = [['1', 'a', '-', '0'], ['2', 'b', '-', '0'], ['3', 'c', '-', '0'], ['4', 'd', '-', '0'],
                 ['5', 'e', '-', '0'], ['6', 'f', '-', '0']]

        parse(lines, tree)

        for line in lines:
            node = tree.nodes[line[0]]
            self.assertEqual(node.name, line[1])
            self.assertEqual(node.rank, ranks[line[3]])

    def test_map_names(self):

        # check if every node is assigned a name
        tree = get_phylo_tree('../../resources/gtdb.tre')
        map_names('../../resources/gtdb.map', tree)
        for n in tree.nodes.values():
            self.assertIsNotNone(n.name)
            self.assertIsNotNone(n.rank)


if __name__ == '__main__':
    unittest.main()
