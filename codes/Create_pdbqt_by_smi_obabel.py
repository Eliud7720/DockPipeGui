import subprocess
import os
import glob

class Conversions():
    def __init__(self):
        self.contator = 0
        self.maxim = 0
        
    def conversions(self, folder, file):
        real_folder = folder + "/"
        os.makedirs(real_folder, exist_ok=True)

        # Create empty lists
        list_smiles = []
        list_ID = []

        # Save ID's and SMiles
        with open(file, 'r') as file:
            lines = file.readlines()

        for line in lines:
            list_smiles.append(line.split()[0])
            list_ID.append(line.split()[1])

        for i, smile in enumerate(list_smiles):
            subprocess.run(["./lib/obabel", "-:" + smile, "-omol2", "-O", real_folder + f'ligand_{i}.mol2', '--gen3d'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            sdf = real_folder + f'ligand_{i}.mol2'
            pdbqt = real_folder + f'ligand_{i}.pdbqt'
            subprocess.run(['./lib/obabel', '-imol2', sdf, '-opdbqt', '-O', pdbqt], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.contator += 1
            os.remove(real_folder + f'ligand_{i}.mol2')
            os.rename(pdbqt, real_folder + list_ID[i] + ".pdbqt")
            yield i
        
    def Maximum(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            self.maxim = len(lines)
