import unittest

from pygan.tree.newick_parser import *


class NewickParserTest(unittest.TestCase):

    def test_read(self):
        with self.assertRaises(FileNotFoundError):
            read('not a file.tre')

    def test_parse(self):

        newick = '((1)2)3'
        tree = parse(newick)

        self.assertListEqual([c.tax_id for c in tree.nodes['1'].children], [])
        self.assertListEqual([c.tax_id for c in tree.nodes['2'].children], ['1'])
        self.assertListEqual([c.tax_id for c in tree.nodes['3'].children], ['2'])

        newick = '(1,2,(3,4)5)6'
        tree = parse(newick)

        self.assertListEqual([c.tax_id for c in tree.nodes['1'].children], [])
        self.assertListEqual([c.tax_id for c in tree.nodes['2'].children], [])
        self.assertListEqual([c.tax_id for c in tree.nodes['3'].children], [])
        self.assertListEqual([c.tax_id for c in tree.nodes['4'].children], [])
        self.assertListEqual([c.tax_id for c in tree.nodes['5'].children], ['3', '4'])
        self.assertListEqual([c.tax_id for c in tree.nodes['6'].children], ['1', '2', '5'])

        newick = '((1,2)3,(4,5)6)7'
        tree = parse(newick)

        self.assertListEqual([c.tax_id for c in tree.nodes['1'].children], [])
        self.assertListEqual([c.tax_id for c in tree.nodes['2'].children], [])
        self.assertListEqual([c.tax_id for c in tree.nodes['3'].children], ['1', '2'])
        self.assertListEqual([c.tax_id for c in tree.nodes['4'].children], [])
        self.assertListEqual([c.tax_id for c in tree.nodes['5'].children], [])
        self.assertListEqual([c.tax_id for c in tree.nodes['6'].children], ['4', '5'])
        self.assertListEqual([c.tax_id for c in tree.nodes['7'].children], ['3', '6'])

    def test_to_newick(self):

        newick = '((1)2)3'
        tree = parse(newick)
        self.assertEqual(to_newick(tree.root), newick)

        newick = '(1,2,(3,4)5)6'
        tree = parse(newick)
        self.assertEqual(to_newick(tree.root), newick)

        newick = '((1,2)3,(4,5)6)7'
        tree = parse(newick)
        self.assertEqual(to_newick(tree.root), newick)

    def test_get_phylo_tree(self):

        newick = read('../../resources/gtdb.tre')
        tree = parse(newick)
        self.assertEqual(to_newick(tree.root), newick)


if __name__ == '__main__':
    unittest.main()
