import sqlite3
import os
from typing import Iterable, Dict


def get_accessions2taxonids(database_path: str, accessions: Iterable[str], key: str = 'Taxonomy') -> Dict[str, int]:
    """
    Connect to megan_map.db and create a dictionary of accessions to taxonomy ids

    :param database_path: path of megan_map.db
    :param accessions: collection of accessions to be mapped
    :param key: What the accession should be mapped to. Taxonomy by default.
    :return: dictionary of accessions to taxonomy ids
    """
    connection = connect(database_path)
    cursor = connection.cursor()
    accessions2taxonids = map_accessions(cursor, accessions, key)
    disconnect(connection)
    return accessions2taxonids


def connect(database_path: str) -> sqlite3.Connection:
    """
    Connect to megan_map.db

    :param database_path: path of megan_map.db
    :return: sqlite3 connection to megan_map.db
    """
    if not os.path.isfile(database_path):
        raise FileNotFoundError('Can not connect to ' + database_path)
    connection = sqlite3.connect(database_path)
    cursor = connection.cursor()
    # raises sqlite3.OperationalError if Accession, Taxonomy or mappings does not exist
    cursor.execute('select Accession, Taxonomy from mappings limit 1')
    return connection


def map_accessions(cursor: sqlite3.Cursor, accessions: Iterable[str], key: str = 'Taxonomy') -> Dict[str, int]:
    """
    Create a dictionary of accessions to taxonomy ids

    :param cursor: sqlite3 cursor to megan_map.db
    :param accessions: collection of accessions to be mapped
    :param key: What the accession should be mapped to. Taxonomy by default.
    :return: dictionary of accessions to taxonomy ids
    """
    cursor.execute(f'select Accession, {key} from mappings where Accession in {tuple(accessions)}')
    return {a: t for a, t in cursor.fetchall()}


def disconnect(connection: sqlite3.Connection):
    """
    Disconnect from megan_map.db

    :param connection: sqlite3 connection to megan_map.db
    """
    connection.close()
