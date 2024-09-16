import subprocess
import os
import glob
from rdkit import Chem

class Conversions():
    def __init__(self):
        self.contator = 0
        self.maxim = 0
        
    def conversions(self, folder, file):
        real_folder = folder + "/"
        os.makedirs(real_folder, exist_ok=True)

        # Smiles list
        smiles_list = []

        supplier = Chem.SDMolSupplier(file)

        for mol in supplier:
            smiles = smiles = Chem.MolToSmiles(mol)
            smiles_list.append(smiles)


        for i, smile in enumerate(smiles_list):
            subprocess.run(["./lib/obabel", "-:" + smile, "-omol2", "-O", real_folder + f'ligand_{i}.mol2', '--gen3d'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            sdf = real_folder + f'ligand_{i}.mol2'
            pdbqt = real_folder + f'ligand_{i}.pdbqt'
            subprocess.run(['./lib/obabel', '-imol2', sdf, '-opdbqt', '-O', pdbqt], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.contator += 1
            os.remove(real_folder + f'ligand_{i}.mol2')
            yield i
        
    def Maximum(self, file):
        list_molecules = []

        supplier = Chem.SDMolSupplier(file)
        for mol in supplier:
            if mol is not None:
                list_molecules.append(mol)

        self.maxim = len(list_molecules)
