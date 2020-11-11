import sys
from typing import List, Optional, Set, Tuple, Union

import pandas as pd
import psycopg2
from sqlalchemy import create_engine

from config import *

dbparams = {
    'dbname': SQL_DBNAME,
    'user': SQL_USER,
    'password': SQL_PASSWORD,
    'host': SQL_HOST,
    'port': SQL_PORT,
}

def initialize_compound_db(ids: List[int],
                           smiles: List[str],
                           spectrum_C: List[List[float]],
                           spectrum_H: List[List[float]],
                           solvents_C: List[str],
                           solvents_H: List[str],
                           mol_blocks: List[str],):
    sql_engine_connect = "postgresql://%s:%s@%s:5432/%s" % (
        dbparams['user'],
        dbparams['password'],
        dbparams['host'],
        dbparams['dbname']
    )

    exp_df = pd.DataFrame({'smiles': smiles, 'spectrum_c': spectrum_C, 'spectrum_h': spectrum_H})
    calculate_df = pd.DataFrame({'id': ids, 'smiles': smiles, 'solvent_c': solvents_C, 'solvent_h': solvents_H, 'romol': mol_blocks})

    engine = create_engine(sql_engine_connect)
    with psycopg2.connect(**dbparams) as conn:
        with conn.cursor() as cur:
            cur.execute("select exists(select * from information_schema.tables where table_name=%s)", ('exp_nmr',))
            if not cur.fetchone()[0]:
                cur.execute("""
                    CREATE TABLE exp_nmr (
                        smiles TEXT PRIMARY KEY,
                        spectrum_c numeric[],
                        spectrum_h numeric[]
                    );
                """)

            cur.execute("select exists(select * from information_schema.tables where table_name=%s)", ('compound',))
            if not cur.fetchone()[0]:
                cur.execute("""
                    CREATE TABLE compound (
                        id INTEGER PRIMARY KEY,
                        smiles TEXT,
                        romol TEXT,
                        status TEXT DEFAULT 'initialized',
                        queued_at TIMESTAMP,
                        node TEXT,
                        coords TEXT,
                        nmr numeric[],
                        solvent_c TEXT,
                        solvent_h TEXT,
                        error TEXT
                    );
                """)

    exp_df.to_sql('exp_nmr',
                  con=engine,
                  index=False,
                  if_exists='append')

    calculate_df.to_sql('compound',
                        con=engine,
                        index=False,
                        if_exists='append')


if __name__ == '__main__':
    nmrshiftdb_pickle = sys.argv[1]
    nmrshiftdb_chiral = pd.read_pickle(nmrshiftdb_pickle)
    smiles = nmrshiftdb_chiral.SMILES.tolist()
    mol_blocks = nmrshiftdb_chiral.ROMol.tolist()
    Cs = nmrshiftdb_chiral.spectrum_C.tolist()
    Hs = nmrshiftdb_chiral.spectrum_H.tolist()
    C_solvents = nmrshiftdb_chiral.Solvent_C.tolist()
    H_solvents = nmrshiftdb_chiral.Solvent_C.tolist()
    ids = nmrshiftdb_chiral['nmrshiftdb2 ID'].tolist()

    initialize_compound_db(ids, smiles, Cs, Hs, C_solvents, H_solvents, mol_blocks)
