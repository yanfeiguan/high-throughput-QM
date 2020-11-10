import os
import tempfile
import subprocess
import logging
import shutil
from datetime import datetime, timedelta

from rdkit import Chem
from rdkit.Chem.rdMolTransforms import GetBondLength

from worker.file_parser import G16Log

logging.getLogger("cclib").setLevel(30)


class QmWorker(object):

    def __init__(self,
                 qm_cmd="g16",
                 wall_time=48,      #2*24*3600
                 nprocs=8,
                 mem='32000mb',
                 scratchdir=os.environ['WORKDIR'],
                 projectdir='/home/yanfeig/nmr/test'):
        """ Class to handle the overall temporary directory management for
        running Gaussian on Eagle """

        self.start_time = datetime.now()
        self.wall_time = timedelta(hours=wall_time)

        self.qm_cmd = qm_cmd
        self.nprocs = nprocs
        self.mem = mem
        self.scratchdir = scratchdir
        self.projectdir = projectdir

    def remaining_wall_time(self, buffer=120):
        buffer = timedelta(seconds=buffer)
        time_now = datetime.now()
        remaining = self.wall_time - (time_now - self.start_time) - buffer
        return remaining.total_seconds()

    def run_qm(self, molblock, solvent='dimethylsulfoxide'):
        """ Given an rdkit.Mol object with an optimized, minimum energy conformer
        ID, write a gaussian input file using openbabel to the scratch folder """

        mol = Chem.MolFromMolBlock(molblock, removeHs=False)
        name = mol.GetProp('_Name')

        with tempfile.TemporaryDirectory(dir=self.scratchdir) as running_dir:
            gjf = os.path.join(running_dir, '{}.gjf'.format(name))
            with tempfile.NamedTemporaryFile('wt', suffix='.sdf', dir=running_dir) as sdf_file:
                writer = Chem.SDWriter(sdf_file.name)
                writer.write(mol)
                writer.close()

                header1 = [
                    #'%MEM={}'.format(self.mem),
                    #'%nprocshared={}'.format(self.nprocs),
                    '#mPW1PW91/6-31+G(d,p) scf=(xqc,maxconventionalcycles=400) nmr=giao nosymm SCRF=(solvent={})'.format(solvent)]

                subprocess.run(
                    ['obabel', sdf_file.name, '-O', gjf, '-xk', '\n'.join(header1)])

            scratch_log = os.path.join(running_dir, '{}.log'.format(name))
            cmd = self.qm_cmd + " < {0} > {1}".format(gjf, scratch_log)

            env = os.environ.copy()
            env['GAUSS_SCRDIR'] = running_dir
            subprocess.run(cmd, shell=True, env=env, timeout=self.remaining_wall_time())
            log = G16Log(scratch_log)
            if log.termination:
                coords = log.coords
                mol = self.assign_coords_to_mol(mol, coords)
                mol_block = Chem.MolToMolBlock(mol)
                nmr = log.NMR
                scf = log.scf / 27.2114
            else:
                raise RuntimeError('Gaussian calculation for NMR failed')

            project_gjf = os.path.join(self.projectdir, '{}.nmr.gjf'.format(name))
            project_log = os.path.join(self.projectdir, '{}.nmr.log'.format(name))

            shutil.move(gjf, project_gjf)
            shutil.move(scratch_log, project_log)
            subprocess.run(['gzip', '-f', project_log, project_gjf])

        return mol_block, nmr, scf

    @staticmethod
    def assign_coords_to_mol(mol, coords):

        conf = mol.GetConformer()

        # Correct the 3D positions of the atoms using the optimized geometry
        for i in range(conf.GetNumAtoms()):
            conf.SetAtomPosition(i, coords[i])
        
        covalent_radii = {'H': .31, 'C': .76, 'N': .71,
                          'O': .66, 'P': 1.07, 'S': 1.05,
                          'F': .57, 'Cl': 1.02, 'Br': 1.20}

        # Check bond lengths
        for bond in mol.GetBonds():
            length = GetBondLength(
                mol.GetConformer(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            
            max_length = (covalent_radii[bond.GetBeginAtom().GetSymbol()] + 
                          covalent_radii[bond.GetEndAtom().GetSymbol()] + 0.4)
            
            assert length <= max_length, "bond greater than maximum covalent length"

        return mol


if __name__ == '__main__':
    worker = QmWorker()
    mol = Chem.MolFromSmiles('C')
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    mol = Chem.AddHs(mol, addCoords=True)
    mol.SetProp('_Name', 'test')
    molblock = Chem.MolToMolBlock(mol)
    mol, nmr, scf = worker.run_qm(molblock)
    print(scf)
