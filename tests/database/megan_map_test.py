import unittest
import sqlite3
import os

import database.megan_map as m


connection = sqlite3.connect('../../megan-map-Jan2021.db')
cursor = connection.cursor()

cursor.execute('select Accession,Taxonomy from mappings where Taxonomy is not null limit 2000')
accessions2taxonids = {a: t for a, t in cursor.fetchall()}

cursor.execute('select Accession, GTDB from mappings where GTDB is not null limit 20')
accessions2gtdbids = {a: g for a, g in cursor.fetchall()}


class MeganMapTests(unittest.TestCase):

    def test_connect(self):
        with self.assertRaises(FileNotFoundError):
            m.connect('not-megan-map.db')
        con = sqlite3.connect('any-other-map.db')
        cur = con.cursor()
        cur.execute("""create table if not exists mappings (
            Accession text primary key,
            nottaxonomy text not null        
        );""")
        con.close()
        with self.assertRaises(sqlite3.OperationalError):
            m.connect('any-other-map.db')
        os.remove('any-other-map.db')

    def test_map_accessions(self):
        con = sqlite3.connect('../../megan-map-Jan2021.db')
        cur = con.cursor()
        accs = accessions2taxonids.keys()
        self.assertDictEqual(m.map_accessions(cur, accs), accessions2taxonids)
        accs = accessions2gtdbids.keys()
        self.assertDictEqual(m.map_accessions(cur, accs, key='GTDB'), accessions2gtdbids)
        with self.assertRaises(sqlite3.OperationalError):
            m.map_accessions(cur, accs, key='not a key')
        con.close()

    def test_get_accessions2taxonids(self):
        con = sqlite3.connect('../../megan-map-Jan2021.db')
        cur = con.cursor()
        accs = accessions2taxonids.keys()
        self.assertDictEqual(m.get_accessions2taxonids('../../megan-map-Jan2021.db', accs), accessions2taxonids)
        accs = accessions2gtdbids.keys()
        self.assertDictEqual(m.get_accessions2taxonids('../../megan-map-Jan2021.db', accs, key='GTDB'),
                             accessions2gtdbids)
        with self.assertRaises(sqlite3.OperationalError):
            m.get_accessions2taxonids('../../megan-map-Jan2021.db', accs, key='not a key')
        con.close()


if __name__ == '__main__':
    unittest.main()
