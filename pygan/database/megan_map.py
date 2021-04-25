import sqlite3
import os
from typing import Dict, List


def fetch_all_taxonids(database_path: str, all_accessions: List[List[str]], key: str = 'Taxonomy') -> List[List[int]]:
    """
    Connect to megan_map.db and create a list of taxonomy ids per read

    :param database_path: path of megan_map.db
    :param all_accessions: collection of accessions per read to be mapped
    :param key: What the accession should be mapped to. Taxonomy by default.
    :return: list of taxonomy ids per read
    """
    connection = connect(database_path)
    all_taxonids = [fetch_ids(connection, accessions, key) for accessions in all_accessions]
    disconnect(connection)
    return all_taxonids


def get_accessions2taxonids(database_path: str, accessions: List[str], key: str = 'Taxonomy') -> Dict[str, int]:
    """
    Connect to megan_map.db and create a dictionary of accessions to taxonomy ids

    :param database_path: path of megan_map.db
    :param accessions: collection of accessions to be mapped
    :param key: What the accession should be mapped to. Taxonomy by default.
    :return: dictionary of accessions to taxonomy ids
    """
    connection = connect(database_path)
    accessions2taxonids = map_accessions2ids(connection, accessions, key)
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
    # raises sqlite3.OperationalError if Accession, Taxonomy or mappings does not exist
    connection.execute('select Accession, Taxonomy from mappings limit 1')
    return connection


def map_accessions2ids(connection: sqlite3.Connection, accessions: List[str], key: str = 'Taxonomy') -> Dict[str, int]:
    """
    Create a dictionary of accessions to taxonomy ids

    :param connection: sqlite3 connection to megan_map.db
    :param accessions: collection of accessions to be mapped
    :param key: What the accession should be mapped to. Taxonomy by default.
    :return: dictionary of accessions to taxonomy ids
    """
    prep = f'select Accession, {key} from mappings where Accession in (\''
    return {a: t for a, t in connection.execute(prep + '\',\''.join(accessions) + '\')')}


def fetch_ids(connection: sqlite3.Connection, accessions: List[str], key: str = 'Taxonomy') -> List[int]:
    """
    Create a list of taxonomy ids from accessions

    :param connection: sqlite3 connection to megan_map.db
    :param accessions: collection of accessions to be mapped
    :param key: What the accession should be mapped to. Taxonomy by default.
    :return: list of taxonomy ids
    """
    prep = f'select {key} from mappings where Accession in (\''
    return [t[0] for t in connection.execute(prep + '\',\''.join(accessions) + '\')')]


def disconnect(connection: sqlite3.Connection):
    """
    Disconnect from megan_map.db

    :param connection: sqlite3 connection to megan_map.db
    """
    connection.close()
