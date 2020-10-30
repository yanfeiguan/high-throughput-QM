from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


class RdkitConformer(object):

    def __init__(self, max_conformers=1000, min_conformers=100):
        self.max_conformers = max_conformers
        self.min_conformers = min_conformers

    def gen_confs(self, smiles):
        """ Embed a molecule in 3D space, optimizing a number of conformers and
        selecting the most stable
        """

        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.rdmolops.AddHs(mol)

        # Use min < 3^n < max conformers, where n is the number of rotatable bonds
        NumRotatableBonds = AllChem.CalcNumRotatableBonds(mol)
        NumConformers = np.clip(3 ** NumRotatableBonds, self.min_conformers,
                                self.max_conformers)

        conformers = AllChem.EmbedMultipleConfs(
            mol, numConfs=int(NumConformers), pruneRmsThresh=0.2, randomSeed=1,
            useExpTorsionAnglePrefs=True, useBasicKnowledge=True)

        def optimize_conformer(conformer):
            prop = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
            ff = AllChem.MMFFGetMoleculeForceField(mol, prop, confId=conformer)
            ff.Minimize()
            return float(ff.CalcEnergy())

        assert conformers, "Conformer embedding failed"

        conformer_energies = np.array([optimize_conformer(conformer) for conformer in conformers])

        return mol, conformer_energies


if __name__ == '__main__':
    ConfGen = RdkitConformer()
    mol, energies = ConfGen.gen_confs('CCCCC')
    print('AAAA')