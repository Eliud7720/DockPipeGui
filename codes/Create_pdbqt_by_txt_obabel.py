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
        print(os.getcwd())
        subprocess.run(['./lib/obabel', '-ismi', file, '-omol2', '-O', real_folder + 'ligand_.mol2', '-m', '--gen3d'], check=True)
        
        mol2_files = glob.glob(real_folder + '*.mol2')
        
        for mol2_file in mol2_files:
            pdbqt_file = os.path.join(folder, os.path.basename(mol2_file).replace('.mol2', '.pdbqt'))
            subprocess.run(['./lib/obabel', '-imol2', mol2_file, '-opdbqt', '-O', pdbqt_file], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.contator += 1
            yield pdbqt_file
        
        # Eliminar archivos .mol2 después de la conversión
        for mol2_file in mol2_files:
            os.remove(mol2_file)
        
    def Maximum(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            self.maxim = len(lines)
