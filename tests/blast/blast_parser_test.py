import unittest
import os

from pygan.blast.blast_parser import *


class BlastParserTest(unittest.TestCase):

    def test_filter_by_top_score(self):

        a0 = ('a0', 0)
        a1 = ('a1', 1)
        a2 = ('a2', 2)
        a3 = ('a3', 3)
        a4 = ('a4', 4)
        a5 = ('a5', 5)
        a6 = ('a6', 6)
        a7 = ('a7', 7)
        a8 = ('a8', 8)
        a9 = ('a9', 9)
        a10 = ('a10', 10)

        read = [a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10]

        accs0 = filter_by_top_score(read, 0)
        self.assertListEqual(accs0, ['a10'])

        accs1 = filter_by_top_score(read, 1)
        self.assertListEqual(accs1, [s for s, f in read])

        accs05 = filter_by_top_score(read, 0.5)
        self.assertListEqual(accs05, [s for s, f in [a5, a6, a7, a8, a9, a10]])

        accs01 = filter_by_top_score(read, 0.1)
        self.assertListEqual(accs01, [s for s, f in [a9, a10]])

    def test_parse_tab(self):

        with open('test.bt', 'w') as f:
            text = ''
            for c in 'abcd':
                for i in range(11):
                    text += 'id_' + c + '\t' + c + str(i) + '.1\t' + str(i) + '\n'
            f.write(text)

        ids0 = [[c + str(i) for i in [10]] for c in 'abcd']
        ids1 = [[c + str(i) for i in range(11)] for c in 'abcd']
        ids05 = [[c + str(i) for i in range(5, 11)] for c in 'abcd']
        ids01 = [[c + str(i) for i in range(9, 11)] for c in 'abcd']

        with open('test.bt', 'r') as f:
            self.assertListEqual(parse_tab(f, 0), ids0)
        with open('test.bt', 'r') as f:
            self.assertListEqual(parse_tab(f, 1), ids1)
        with open('test.bt', 'r') as f:
            self.assertListEqual(parse_tab(f, 0.5), ids05)
        with open('test.bt', 'r') as f:
            self.assertListEqual(parse_tab(f, 0.1), ids01)

        with open('test.bt', 'r') as f:
            parsed03 = parse_tab(f, 0.3)
        with open('test.bt', 'r') as f:
            parsed07 = parse_tab(f, 0.7)

        self.assertEqual(len(parsed03), len(parsed07))
        self.assertNotEqual(len(parsed03[0]), len(parsed07[0]))

        os.remove('test.bt')


if __name__ == '__main__':
    unittest.main()
