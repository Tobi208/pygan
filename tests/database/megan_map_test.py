import unittest
import sqlite3

import database.megan_map as m


connection = sqlite3.connect('../../megan-map-Jan2021.db')
cursor = connection.cursor()
cursor.execute('select Accession from mappings limit 2000')
accessions = [x[0] for x in cursor.fetchall()]
accessions2taxonids = {}
for accession in accessions:
    cursor.execute(f'select Taxonomy from mappings where Accession=\'{accession}\'')
    accessions2taxonids[accession] = cursor.fetchone()[0]


class MeganMapTests(unittest.TestCase):

    def test_map_accessions(self):
        con = sqlite3.connect('../../megan-map-Jan2021.db')
        cur = con.cursor()
        accs = accessions2taxonids.keys()
        self.assertDictEqual(m.map_accessions(cur, accs), accessions2taxonids)
        con.close()

    def test_get_accessions2taxonids(self):
        con = sqlite3.connect('../../megan-map-Jan2021.db')
        cur = con.cursor()
        accs = accessions2taxonids.keys()
        self.assertDictEqual(m.get_accessions2taxonids('../../megan-map-Jan2021.db', accs), accessions2taxonids)
        con.close()


if __name__ == '__main__':
    unittest.main()
