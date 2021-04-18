import unittest
import sqlite3
import os

import pygan.database.megan_map as m


connection = sqlite3.connect('../../resources/megan-map-Jan2021.db')
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
        con = sqlite3.connect('../../resources/megan-map-Jan2021.db')
        cur = con.cursor()
        accs = list(accessions2taxonids.keys())
        self.assertDictEqual(m.map_accessions2ids(cur, accs), accessions2taxonids)
        accs = list(accessions2gtdbids.keys())
        self.assertDictEqual(m.map_accessions2ids(cur, accs, key='GTDB'), accessions2gtdbids)
        with self.assertRaises(sqlite3.OperationalError):
            m.map_accessions2ids(cur, accs, key='not a key')
        con.close()

    def test_get_accessions2taxonids(self):
        accs = list(accessions2taxonids.keys())
        self.assertDictEqual(m.get_accessions2taxonids('../../resources/megan-map-Jan2021.db', accs),
                             accessions2taxonids)
        accs = list(accessions2gtdbids.keys())
        self.assertDictEqual(m.get_accessions2taxonids('../../resources/megan-map-Jan2021.db', accs, key='GTDB'),
                             accessions2gtdbids)
        with self.assertRaises(sqlite3.OperationalError):
            m.get_accessions2taxonids('../../resources/megan-map-Jan2021.db', accs, key='not a key')

    def test_real_accessions(self):

        accessions = ['EXZ85837']
        a2id = m.get_accessions2taxonids('../../resources/megan-map-Jan2021.db', accessions)
        self.assertDictEqual(a2id, {'EXZ85837': 1339273})

        accessions = ['EXZ85837', 'WP_065538752', 'MBE5744624']
        a2id = m.get_accessions2taxonids('../../resources/megan-map-Jan2021.db', accessions)
        self.assertDictEqual(a2id, {
            'EXZ85837': 1339273,
            'WP_065538752': 1796613
        })

    def test_fetch_ids(self):

        con = sqlite3.connect('../../resources/megan-map-Jan2021.db')
        cur = con.cursor()

        accessions = ['EXZ85837']
        ids = m.fetch_ids(cur, accessions)
        self.assertListEqual(ids, [1339273])

        accessions = ['EXZ85837', 'WP_065538752']
        ids = m.fetch_ids(cur, accessions)
        self.assertListEqual(ids, [1339273, 1796613])

        accessions = ['EXZ85837', 'WP_065538752', 'MBE5744624']
        ids = m.fetch_ids(cur, accessions)
        self.assertListEqual(ids, [1339273, 1796613])


if __name__ == '__main__':
    unittest.main()
