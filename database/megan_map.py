import sqlite3
from typing import List


class MeganMap:
    """
    Interface between megan map database and the program
    """

    def __init__(self):
        self.connection = None
        self.cursor = None
        self.accension2taxonomy = {}

    def connect(self, db_file: str):
        """
        Connect megan map to database file

        :param db_file: database file containing megan map
        """
        self.connection = sqlite3.connect(db_file)
        self.cursor = self.connection.cursor()

    def load_accession(self, accession: str) -> None:
        """
        Select taxonomy from database by accession and add it to dictionary

        :param accession: accession mapped to taxonomy
        """
        self.cursor.execute(f'select Taxonomy from mappings where Accession=\'{accession}\'')
        self.accension2taxonomy[accession] = self.cursor.fetchone()[0]

    def load_accessions(self, accessions: List[str]) -> None:
        """
        Select taxonomies from database by accessions and add them to dictionary

        :param accessions: list of accessions mapped to taxonomies
        :return:
        """
        for accession in accessions:
            self.load_accession(accession)

    def get_taxonomy(self, accension: str) -> str:
        """
        Retrieve taxonomy by accession from dictionary or load from database

        :param accension: accession mapped to taxonomy
        :return: taxonomy mapped to accession
        """
        if accension not in self.accension2taxonomy:
            self.load_accession(accension)
        return self.accension2taxonomy[accension]

    def get_taxonids(self, accensions: List[str]) -> List[str]:
        """
        Retrieve taxonomies by accessions from dictionary or load from database

        :param accensions: accessions mapped to taxonomies
        :return: taxonomies mapped to accessions
        """
        return list(map(self.get_taxonomy, accensions))
