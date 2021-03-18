import unittest
import os
from tree.map_parser import *
from tree.newick_parser import parse as nwk_parse, get_phylo_tree


class MapParserTest(unittest.TestCase):

    def test_readlines(self):

        with self.assertRaises(FileNotFoundError):
            readlines('not a file.map')

        with open('test.map', 'w') as f:
            f.write('1\ta\tblabla\n2\tb\tblabla\n3\tc\tblabla\n')
        self.assertListEqual(readlines('test.map'), [['1', 'a'], ['2', 'b'], ['3', 'c']])
        os.remove('test.map')

    def test_parse(self):

        tree = nwk_parse('(1,2,(3,4)5)6')
        lines = [['1', 'a'], ['2', 'b'], ['3', 'c'], ['4', 'd'], ['5', 'e'], ['6', 'f']]

        parse(lines, tree)

        for line in lines:
            self.assertEqual(tree.nodes[line[0]].name, line[1])

    def test_map_names(self):

        # check if every node is assigned a name
        tree = get_phylo_tree('../../resources/gtdb.tre')
        map_names('../../resources/gtdb.map', tree)
        for k, n in tree.nodes.items():
            self.assertIsNotNone(n.name)


if __name__ == '__main__':
    unittest.main()
