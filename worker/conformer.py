from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


class RdkitConformer(object):

    def __init__(self, max_conformers=20, min_conformers=10, e_threshold=2.5, rmsd_threshold=0.4):
        self.max_conformers = max_conformers
        self.min_conformers = min_conformers
        self.e_threshold = e_threshold
        self.rmsd_threshold = rmsd_threshold

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

        conformer_energies = [optimize_conformer(conformer) for conformer in conformers]
        conformer_ids = range(len(conformer_energies))

        # filter through energy threshold and rmsd
        conformer_energies, conformer_ids = zip(*sorted(zip(conformer_energies, conformer_ids)))
        n_conformers_valid = len([x for x in conformer_energies if x - conformer_energies[0] <= self.e_threshold])
        conformer_ids = conformer_ids[:n_conformers_valid]
        conformer_ids = self.postrmsd(mol, conformer_ids)

        return mol, conformer_ids

    # filter conformers based on geometric RMS
    def postrmsd(self, mol, conformer_ids):

        conformer_ids_baned = []
        for i, confId_i in enumerate(conformer_ids):
            if confId_i in conformer_ids_baned:
                continue
            for confId_j in conformer_ids[i+1:]:
                if confId_j in conformer_ids_baned:
                    continue
                rmsd = AllChem.GetBestRMS(mol, mol, prbId=confId_i, refId=confId_j)
                if rmsd < self.rmsd_threshold:
                    conformer_ids_baned.append(confId_j)

        return [x for x in conformer_ids if x not in conformer_ids_baned]


if __name__ == '__main__':
    ConfGen = RdkitConformer()
    mol, energies = ConfGen.gen_confs('CCCCC')
    mol.SetProp('_Name', 'test')