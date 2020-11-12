import time
import socket
import subprocess

from sqlalchemy import create_engine
import psycopg2
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd

from worker import RdkitConformer, XtbOptimizer, QmWorker
from config import *

dbparams = {
    'dbname': SQL_DBNAME,
    'user': SQL_USER,
    'password': SQL_PASSWORD,
    'host': SQL_HOST,
    'port': SQL_PORT,
}

R = 0.001987
T = 298.15

opt_columns = ['mol_opt', 'e', 'h', 'g', 'error_opt', 'status_opt']
qm_columns = ['mol_qm', 'nmr', 'scf', 'error_qm', 'status_qm']
#from bde.gaussian import GaussianRunner


def boltzmann(energies):
    energies_e = np.exp(-energies/(R*T))
    dist = energies_e/np.nansum(energies_e)
    return dist


class JobWrapper(object):

    def __init__(self,
                 conformer_generator=RdkitConformer(),
                 optimizer=XtbOptimizer(),
                 qm_worker=QmWorker(),
                 optimization_e_threshold=3,
                 optimization_rmsd=0.4):

        self.conformer_generator = conformer_generator
        self.optimizer = optimizer
        self.qm_worker = qm_worker

        self.optimization_e_threshold = optimization_e_threshold
        self.optimization_rmsd = optimization_rmsd
        self.sql_engine_connect = "postgresql://%s:%s@%s:5432/%s" % (
            dbparams['user'],
            dbparams['password'],
            dbparams['host'],
            dbparams['dbname']
        )
        self.init_conformers()
        self.conn = psycopg2.connect(**dbparams)

    @staticmethod
    def init_conformers():
        with psycopg2.connect(**dbparams) as conn:
            with conn.cursor() as cur:
                cur.execute("select exists(select * from information_schema.tables where table_name=%s)", ('conformers',))
                if not cur.fetchone()[0]:
                    cur.execute("""
                        CREATE table conformers (
                            cid integer, 
                            confid integer,
                            mol_opt text,
                            e numeric,
                            h numeric,
                            g numeric,
                            error_opt text,
                            status_opt text default 'initialized',
                            mol_qm text,
                            nmr numeric[],
                            scf numeric,
                            error_qm text,
                            status_qm text default 'initialized',
                            primary key (cid, confid)
                        )
                    """)

    def batch_optimization(self, conformers, conformer_ids):
        optimized_mols = []
        for finished, (confId, mol_block) in enumerate(zip(conformer_ids, conformers)):
            try:
                mol_opt_block, E, H, G = self.optimizer.run(mol_block)
                optimized_mols.append([confId, mol_opt_block, E, H, G, np.nan, 'finished'])

            except subprocess.TimeoutExpired:
                break
            except Exception as e:
                optimized_mols.append([confId, mol_block, np.nan, np.nan, np.nan, str(e), 'failed'])

        if finished + 1 < len(conformers):      # timeout
            optimized_mols.extend([[conformer_ids[i], conformers[i], np.nan, np.nan, np.nan, np.nan, 'initialized'] for i in range(finished, len(conformer_ids))])

        optimized_mols = pd.DataFrame(optimized_mols, columns=['confid', 'mol_opt', 'e', 'h', 'g', 'error_opt', 'status_opt']).astype({'e': float, 'h': float, 'h': float})

        return optimized_mols, finished + 1 < len(conformers)

    def postprocess(self, optimized_mols):
        # FIXME unit
        optimized_mols = optimized_mols.copy()
        min_E, min_H, min_G = optimized_mols[optimized_mols.status_qm == 'finished'][['e', 'h', 'g']].min().values
        for e, min_e in zip(['e', 'h', 'g'], [min_E, min_H, min_G]):
            optimized_mols['min_' + e] = optimized_mols[e].apply(lambda x: (x - min_e) * 23.0609)

        optimized_mols = optimized_mols.sort_values(by='min_g')
        coords = optimized_mols[(optimized_mols.status_qm == 'finished') & (optimized_mols.g == min_G)].iloc[0]['mol_opt']
        optimized_mols['mol_opt'] = optimized_mols.mol_opt.apply(lambda x: Chem.MolFromMolBlock(x))
        # postrmsd filter
        baned = np.zeros(len(optimized_mols)).astype(bool)
        for i in range(len(optimized_mols)):
            if not pd.isna(optimized_mols.iloc[i]['error_opt']) or baned[i]:
                continue
            for j in range(i+1, len(optimized_mols)):
                if not pd.isna(optimized_mols.iloc[j]['error_opt']) or baned[j]:
                    continue
                rmsd = AllChem.GetBestRMS(optimized_mols.iloc[i]['mol_opt'], optimized_mols.iloc[j]['mol_opt'])
                if rmsd < self.optimization_rmsd:
                    baned[j] = True
        coords_checked_df = optimized_mols[~baned]

        nmr = coords_checked_df.nmr.apply(np.array).tolist()
        boltzmann_weight = boltzmann(coords_checked_df['min_g'].values)

        weighted_nmr = sum([w*n for w, n in zip(boltzmann_weight, nmr) if not np.any(np.isnan(n))]).tolist()

        return weighted_nmr, coords

    def batch_qm(self, mol_blocks, conformer_ids):
        qm_mols = []
        for finished, (mol_block, confId) in enumerate(zip(mol_blocks, conformer_ids)):
            try:
                mol_block, nmr, scf = self.qm_worker.run_qm(mol_block, level_of_theory=LEVEL_OF_THEORY,
                                                            additional_params=ADDITIONAL_PARAMS, solvent=SOLVENT)
                qm_mols.append([confId, mol_block, nmr, scf, np.nan, 'finished'])
            except subprocess.TimeoutExpired:
                break
            except Exception as e:
                qm_mols.append([confId, mol_block, np.nan, np.nan, str(e), 'failed'])

        if finished + 1 < len(mol_blocks):      # timeout
            qm_mols.extend([[conformer_ids[i], mol_blocks[i], np.nan, np.nan, np.nan, 'initialized'] for i in range(finished, len(conformer_ids))])

        qm_mols = pd.DataFrame(qm_mols, columns=['confid', 'mol_qm', 'nmr', 'scf', 'error_qm', 'status_qm'])
        return qm_mols, finished + 1 < len(mol_blocks)

    def update_conformers(self, conformers_df, cid):
        conformers_df['cid'] = cid
        with psycopg2.connect(**dbparams) as conn:
            with conn.cursor() as cur:
                cur.execute("""
                DELETE FROM conformers
                where conformers.cid = %s
                """, (cid, ))

        engine = create_engine(self.sql_engine_connect)
        conformers_df.to_sql(
            'conformers',
            con=engine,
            index=False,
            if_exists='append'
        )

    @staticmethod
    def update_one_entry(cid, status, nmr=None, coords=None, error=None):
        with psycopg2.connect(**dbparams) as conn:
            with conn.cursor() as cur:
                cur.execute("""
                    UPDATE compound
                    SET status = %s, coords = %s, nmr = %s, error = %s 
                    WHERE id = %s;""", (status, coords, nmr, error, cid,)
                )

    @staticmethod
    def grab_one_entry():
        with psycopg2.connect(**dbparams) as conn:
            with conn.cursor() as cur:
                cur.execute("""
                       WITH cte AS (
                       SELECT id, smiles, status
                       FROM compound
                       WHERE status = 'initialized' OR status = 'timeout'
                       ORDER BY id
                       LIMIT 1
                       FOR UPDATE
                       )

                       UPDATE compound SET status = 'running',
                                           queued_at = CURRENT_TIMESTAMP,
                                           node = %s
                       FROM cte
                       WHERE compound.id = cte.id
                       RETURNING cte.id, cte.smiles, cte.status;
                       """, (socket.gethostname(),))

                cid, smiles, status = cur.fetchone()

        conn.close()

        return cid, smiles, status

    @staticmethod
    def grab_conformers(cid):

        def float_nmr(nmr):
            return [float(x) for x in nmr]

        with psycopg2.connect(**dbparams) as conn:
            sql = "SELECT * FROM conformers where cid = {}".format(cid)

            conformer_df = pd.read_sql_query(sql, conn)

        conformer_df['nmr'] = conformer_df.nmr.apply(lambda x: float_nmr(x) if x is not None else x)
        return conformer_df

    def run(self):

        cid, smiles, status = self.grab_one_entry()
        print(cid)

        if status == 'initialized':
            mol, conformer_ids = self.conformer_generator.gen_confs(smiles)
            print(conformer_ids)
            mol.SetProp('_Name', str(cid))

            conformers = []
            for confId in conformer_ids:
                mol.SetProp('_Name', str(cid) + '_conf{}'.format(confId))
                conformers.append(Chem.MolToMolBlock(mol, confId=confId))

            done_opt_df, timeout = self.batch_optimization(conformers, conformer_ids)

            if timeout:
                self.update_one_entry(cid, 'timeout')
                self.update_conformers(done_opt_df, cid)
                return cid, False

            if 'finished' not in done_opt_df.status_opt.tolist():
                self.update_one_entry(cid, 'failed', error='no optimized conformers found')
                return cid, False

            optimized_mols, conformer_ids = zip(*done_opt_df[done_opt_df.status_opt == 'finished'][['mol_opt', 'confid']].values.tolist())

        elif status == 'timeout':
            conformer_df = self.grab_conformers(cid)

            optimized_mols_df = conformer_df[conformer_df.status_opt != 'initialized']
            if len(optimized_mols_df) != len(conformer_df):
                conformers, conformer_ids = zip(*conformer_df[conformer_df.status_opt == 'initialized'][['mol_opt', 'confid']].values.tolist())
                done_opt_df, timeout = self.batch_optimization(conformers, conformer_ids)
                done_opt_df = pd.concat([optimized_mols_df, done_opt_df])
                if timeout:
                    self.update_one_entry(cid, 'timeout')
                    self.update_conformers(done_opt_df, cid)
                    return cid, False
                done_opt_df.status_qm = 'initialized'
            else:
                done_opt_df = optimized_mols_df

            done_qm_df = done_opt_df[done_opt_df.status_qm != 'initialized'].drop(opt_columns, axis=1)
            optimized_mols, conformer_ids = zip(*done_opt_df[(done_opt_df.status_qm == 'initialized') & (done_opt_df.status_opt == 'finished')][['mol_opt', 'confid']].values.tolist())
            done_opt_df = done_opt_df.drop(qm_columns + ['cid'], axis=1)

        qm_mols_df, timeout = self.batch_qm(optimized_mols, conformer_ids)

        if status == 'timeout' and len(done_qm_df) != 0:
            done_qm_df = pd.concat([done_qm_df, qm_mols_df], ignore_index=True)
        else:
            done_qm_df = qm_mols_df

        done_df = done_opt_df.join(done_qm_df.set_index('confid'), on='confid')
        self.update_conformers(done_df, cid)

        if timeout:
            self.updated_one_entry(cid, 'timeout')
            return cid, False

        if 'finished' not in done_df.status_qm.tolist():
            self.update_one_entry(cid, 'failed', error='no qm calculations normally terminated')
            return cid, False

        weighted_nmr, coords = self.postprocess(done_df)
        self.update_one_entry(cid, 'finished', weighted_nmr, coords)

        return cid, True
                

if __name__ == "__main__":

    conformer_generator = RdkitConformer(max_conformers=MAX_CONFORMERS,
                                         min_conformers=MIN_CONFORMERS,
                                         e_threshold=CONFORMER_E_THRESHOLD,
                                         rmsd_threshold=CONFORMER_RMSD_THRESHOLD)

    optimizer = XtbOptimizer(xtb_path=XTB_CMD,
                             scratchdir=SCRATCHDIR,
                             projectdir=PROJECTDIR,
                             wall_time=WALLTIME)

    qm_worker = QmWorker(qm_cmd=QM_CMD,
                         wall_time=WALLTIME,
                         nprocs=NPROCS,
                         mem=MEM,
                         scratchdir=SCRATCHDIR,
                         projectdir=PROJECTDIR)

    job_wrapper = JobWrapper(conformer_generator=conformer_generator,
                             optimizer=optimizer,
                             qm_worker=qm_worker)

    import random

    while True:  # Time in days
        try:
            cid, finished = job_wrapper.run()
            print(cid, finished)
        except psycopg2.OperationalError:
            time.sleep(5 + random.uniform(0, 60))