import psycopg2
import time
import socket

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd

from worker import RdkitConformer, XtbOptimizer, QmWorker

dbparams = {
    'dbname': 'nmrtest',
    'user': 'nmruser',
    'password': '123456',
    'host': 'localhost',
    'port': 5432,
}


#from bde.gaussian import GaussianRunner

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

    def batch_optimization(self, mol, conformer_ids):
        optimized_mols = []
        cid = mol.GetProp('_Name')
        for confId in conformer_ids:
            mol.SetProp('_Name', cid + '_conf{}'.format(confId))
            mol_block = Chem.MolToMolBlock(mol, confId=confId)
            #try:
            mol_opt, E, H, G = self.optimizer.run(mol_block)
            optimized_mols.append([confId, mol_opt, E, H, G, np.nan])
            #except Exception as e:
                #optimized_mols.append([confId, np.nan, np.nan, np.nan, np.nan, e])

        optimized_mols = pd.DataFrame(optimized_mols, columns=['confId', 'mol', 'E', 'H', 'G', 'Error']).astype({'E': float, 'H': float, 'G': float})

        optimized_failed = optimized_mols[~pd.isna(optimized_mols.Error)]
        optimized_success = optimized_mols[pd.isna(optimized_mols.Error)]
        print("Optimization for {} {}: failed: {}, success {}".format(Chem.MolToSmiles(mol), cid, len(optimized_failed), len(optimized_success)))

        # FIXME unit
        min_E, min_H, min_G = optimized_mols[['E', 'H', 'G']].min().values
        for e, min_e in zip(['E', 'H', 'G'], [min_E, min_H, min_G]):
            optimized_success[e] = optimized_success[e].apply(lambda x: (x - min_e) * 23.0609)

        # postrmsd filter
        baned = np.zeros(len(optimized_success)).astype(bool)
        for i in range(len(optimized_success)):
            if baned[i]:
                continue
            for j in range(i+1, len(optimized_success)):
                if baned[j]:
                    continue
                rmsd = AllChem.GetBestRMS(optimized_success.iloc[i]['mol'], optimized_success.iloc[j]['mol'])
                if rmsd < self.optimization_rmsd:
                    baned[j] = True
        optimized_success = optimized_success[~baned]

        optimized_mols['mol'] = optimized_mols.mol.apply(lambda x: Chem.MolToMolBlock(x))
        self.optimizer.save_summary(optimized_mols, cid)

        return optimized_success

    def run(self):

        with psycopg2.connect(**dbparams) as conn:
            with conn.cursor() as cur:
                cur.execute("""
                WITH cte AS (
                SELECT id, smiles
                FROM compound
                WHERE status = 'initialized'
                ORDER BY id
                LIMIT 1
                FOR UPDATE
                )
    
                UPDATE compound SET status = 'running',
                                    queued_at = CURRENT_TIMESTAMP,
                                    node = %s
                FROM cte
                WHERE compound.id = cte.id
                RETURNING compound.id, compound.smiles;
                """, (socket.gethostname(),))

                cid, smiles = cur.fetchone()

        conn.close()

        #try:
        mol, conformer_ids = self.conformer_generator.gen_confs(smiles)
        mol.SetProp('_Name', str(cid))
        optimized_success = self.batch_optimization(mol, conformer_ids)

        runner = GaussianRunner(smiles, cid)
        molstr, enthalpy, freeenergy, scfenergy, log = runner.process()

        with psycopg2.connect(**dbparams) as conn:
            with conn.cursor() as cur:
                cur.execute("""
                UPDATE compound
                SET status = 'finished', 
                    mol = %s, enthalpy = %s, 
                    freeenergy = %s, scfenergy= %s,
                    run_at = CURRENT_TIMESTAMP,
                    logfile = %s
                WHERE id = %s;""",
                (molstr, enthalpy, freeenergy, scfenergy, log, cid))

        conn.close()

        '''
        except Exception as ex:
            with psycopg2.connect(**dbparams) as conn:
                with conn.cursor() as cur:
                    cur.execute("""
                    UPDATE compound
                    SET status = 'error', 
                        error = %s, 
                        run_at = CURRENT_TIMESTAMP
                    WHERE id = %s;""", (str(ex), cid))

            conn.close()
        '''

        return cid
                

if __name__ == "__main__":

    conformer_generator = RdkitConformer()
    optimizer = XtbOptimizer(xtb_path='xtb', scratchdir='/tmp/yanfei/', projectdir='/Users/yanfei/Projects/NMR/workflow/test')
    qm_worker = QmWorker()

    job_wrapper = JobWrapper(conformer_generator=conformer_generator,
                             optimizer=optimizer,
                             qm_worker=qm_worker)
    
    start_time = time.time()
    job_wrapper.run()

    '''
    # Add a random delay to avoid race conditions at the start of the job
    time.sleep(random.uniform(0, 1*60))
    
    while (time.time() - start_time) < (86400 * 9):  # Time in days
        

        try:
            run_optimization()
            
        except psycopg2.OperationalError:
            time.sleep(5 + random.uniform(0, 60))
    '''