import sqlite3
from typing import Iterable, Dict


def get_accessions2taxonids(database_path: str, accessions: Iterable[str]) -> Dict[str, int]:
    """
    Connect to megan_map.db and create a dictionary of accessions to taxonomy ids

    :param database_path: path of megan_map.db
    :param accessions: collection of accessions to be mapped
    :return: dictionary of accessions to taxonomy ids
    """
    connection = connect(database_path)
    cursor = connection.cursor()
    accessions2taxonids = map_accessions(cursor, accessions)
    disconnect(connection)
    return accessions2taxonids


def connect(database_path: str) -> sqlite3.Connection:
    """
    Connect to megan_map.db

    :param database_path: path of megan_map.db
    :return: sqlite3 connection to megan_map.db
    """
    return sqlite3.connect(database_path)


def map_accessions(cursor: sqlite3.Cursor, accessions: Iterable[str]) -> Dict[str, int]:
    """
    Create a dictionary of accessions to taxonomy ids

    :param cursor: sqlite3 cursor to megan_map.db
    :param accessions: collection of accessions to be mapped
    :return: dictionary of accessions to taxonomy ids
    """
    return {accession: inquire(cursor, accession) for accession in accessions}


def inquire(cursor: sqlite3.Cursor, accession: str) -> int:
    """
    Send a query to megan_map.db

    :param cursor: sqlite3 cursor to megan_map.db
    :param accession: accession to be queried
    :return: corresponding taxonomy id
    """
    cursor.execute(f'select Taxonomy from mappings where Accession=\'{accession}\'')
    return cursor.fetchone()[0]


def disconnect(connection: sqlite3.Connection):
    """
    Disconnect from megan_map.db

    :param connection: sqlite3 connection to megan_map.db
    """
    connection.close()
