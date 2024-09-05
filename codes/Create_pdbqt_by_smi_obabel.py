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

        # Create a temporal file withouth ID's
        temporal_file = real_folder + 'temp.txt'
        with open(temporal_file, 'w') as file:
            for smile in list_smiles:
                file.write(smile + "\n")


        
        subprocess.run(['./lib/obabel', '-ismi', temporal_file, '-omol2', '-O', real_folder + 'ligand_.mol2', '-m', '--gen3d'], check=True)
        

        mol2_files = glob.glob(real_folder + '*.mol2')
        mol2_files = sorted(mol2_files, key=lambda x: int(os.path.basename(x).split('_')[-1].replace('.mol2', '')))
        
        for mol2_file in mol2_files:
            pdbqt_file = os.path.join(folder, os.path.basename(mol2_file).replace('.mol2', '.pdbqt'))
            subprocess.run(['./lib/obabel', '-imol2', mol2_file, '-opdbqt', '-O', pdbqt_file], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.contator += 1
            yield pdbqt_file

        for mol2_file in mol2_files:
            os.remove(mol2_file)
        
        os.remove(temporal_file)

        pdbqt_files = glob.glob(real_folder + '*.pdbqt')
        pdbqt_files = sorted(pdbqt_files, key=lambda x: int(os.path.basename(x).split('_')[-1].replace('.pdbqt', '')))


        for i, pdbqt in enumerate(pdbqt_files):
            os.rename(pdbqt, real_folder + list_ID[i] + ".pdbqt")

        
    def Maximum(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            self.maxim = len(lines)
