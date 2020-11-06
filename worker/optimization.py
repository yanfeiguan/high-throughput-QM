import os
import shutil
import subprocess
import tempfile
from datetime import datetime, timedelta

import numpy as np
from rdkit import Chem

from worker.file_parser import XtbLog


class XtbOptimizer(object):

    def __init__(self,
                 xtb_path='xtb',
                 scratchdir='',
                 projectdir='',
                 wall_time=48,  #2*24*3600
                 ):
        self.start_time = datetime.now()
        self.wall_time = timedelta(hours=wall_time)

        self.xtb_path = xtb_path
        self.scratchdir = os.path.join(scratchdir)
        self.projectdir = os.path.join(projectdir, 'xtb')

        self.cmd = xtb_path

        if not os.path.isdir(self.projectdir):
            os.mkdir(self.projectdir)

    def remaining_wall_time(self, buffer=120):
        buffer = timedelta(seconds=buffer)
        time_now = datetime.now()
        remaining = self.wall_time - (time_now - self.start_time) - buffer
        return remaining.total_seconds()

    def optimize(self,
                 mol_block: str = '') -> Chem.Mol:

        mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        name = mol.GetProp('_Name')

        with tempfile.TemporaryDirectory(dir=self.scratchdir) as running_dir:
            writer = Chem.SDWriter(os.path.join(running_dir, '{}.sdf'.format(name)))
            writer.write(mol)
            writer.close()

            with open(os.path.join(running_dir, '{}.xtbopt.log'.format(name)), 'w') as opt_log:
                subprocess.call([self.cmd,
                                 os.path.join(running_dir, '{}.sdf'.format(name)),
                                 '-opt',
                                 '--ceasefiles',
                                 '--namespace', os.path.join(running_dir, name)],
                                stdout=opt_log, stderr=opt_log,
                                timeout=self.remaining_wall_time())

            opt_sdf = os.path.join(running_dir, '{}.xtbopt.sdf'.format(name))
            new_opt_sdf = os.path.join(self.projectdir, '{}.xtbopt.sdf'.format(name))
            shutil.move(opt_sdf, new_opt_sdf)

            opt_log = os.path.join(running_dir, '{}.xtbopt.log'.format(name))
            new_opt_log = os.path.join(self.projectdir, '{}.xtbopt.log'.format(name))
            shutil.move(opt_log, new_opt_log)

        mol_opt = Chem.SDMolSupplier(new_opt_sdf, removeHs=False)[0]
        mol_opt.SetProp('_Name', name + '_xtbopt')
        subprocess.run(['gzip', new_opt_log, new_opt_sdf])

        return mol_opt

    def valid_freq(self,
                   mol_block: str = ''):

        mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        name = mol.GetProp('_Name')

        with tempfile.TemporaryDirectory(dir=self.scratchdir) as running_dir:
            writer = Chem.SDWriter(os.path.join(running_dir, '{}.sdf'.format(name)))
            writer.write(mol)
            writer.close()

            with open(os.path.join(running_dir, '{}.xtbfreq.log'.format(name)), 'w') as freq_log:
                subprocess.call([self.cmd,
                                 os.path.join(running_dir, '{}.sdf'.format(name)),
                                 '--hess',
                                 '--ceasefiles',
                                 '--namespace', os.path.join(running_dir, name)],
                                stdout=freq_log, stderr=freq_log,
                                timeout=self.remaining_wall_time())

            freq_log = os.path.join(running_dir, '{}.xtbfreq.log'.format(name))
            log = XtbLog(freq_log)

            new_freq_log = os.path.join(self.projectdir, '{}.xtbfreq.log'.format(name))
            shutil.move(freq_log, new_freq_log)
            subprocess.run(['gzip', new_freq_log])

        if log.termination:
            peaks = log.wavenum
            if np.min(peaks) < 0:
                return None, None, None
            else:
                return log.E, log.H, log.G
        else:
            return None, None, None

    def save_summary(self, df, cid):
        summary_csv = os.path.join(self.projectdir, '{}_xtbopt.csv'.format(cid))
        df.to_csv(summary_csv)
        subprocess.run(['gzip', summary_csv])

    def run(self,
            mol_block: str = ''):
        mol_opt = self.optimize(mol_block)
        mol_opt_block = Chem.MolToMolBlock(mol_opt)
        E, H, G = self.valid_freq(mol_opt_block)

        if E:
            return mol_opt_block, E, H, G
        else:
            raise RuntimeError('optimization with GFN_XTB failed')


if __name__ == '__main__':
    m = Chem.MolFromSmiles("CCC")
    m.SetProp('_Name', 'CCC')

    from rdkit.Chem import AllChem

    AllChem.EmbedMolecule(m, AllChem.ETKDG())
    m = Chem.AddHs(m, addCoords=True)

    mol_block = Chem.MolToMolBlock(m)

    optimizer = XtbOptimizer(xtb_path='xtb', scratchdir='/tmp/yanfei/', projectdir='/Users/yanfei/Projects/NMR/workflow/test')
    mol = optimizer.optimize(mol_block)

    mol_block = Chem.MolToMolBlock(mol)
    E, H, G = optimizer.valid_freq(mol_block)

    print('AAA')
