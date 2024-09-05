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

        subprocess.run(['./lib/obabel', '-isdf', file, '-omol2', '-O', real_folder + 'ligand_.mol2', '-m', '--gen3d'], check=True)

        mol2_files = glob.glob(real_folder + '*.mol2')

        # Ordenar archivos por nombre para asegurar orden correcto
        mol2_files = sorted(mol2_files, key=lambda x: int(os.path.basename(x).split('_')[-1].replace('.mol2', '')))

        for mol2_file in mol2_files:
            pdbqt_file = os.path.join(folder, os.path.basename(mol2_file).replace('.mol2', '.pdbqt'))
            subprocess.run(['./lib/obabel', '-imol2', mol2_file, '-opdbqt', '-O', pdbqt_file], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.contator += 1
            yield pdbqt_file

        for mol2_file in mol2_files:
            os.remove(mol2_file)
        
    def Maximum(self, file):
        list_molecules = []

        supplier = Chem.SDMolSupplier(file)
        for mol in supplier:
            if mol is not None:
                list_molecules.append(mol)

        self.maxim = len(list_molecules)
