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

        with open('errors.txt', 'w') as error_file:
            error_file.write("Moléculas que fallaron al generar coordenadas 3D:\n")
        
        result = subprocess.run(['./lib/obabel', '-ismi', file, '-omol2', '-O', real_folder + 'ligand_.mol2', '-m', '--gen3d'],
                                capture_output=True, text=True)
        

        mol2_files = glob.glob(real_folder + '*.mol2')
        mol2_files = sorted(mol2_files, key=lambda x: int(os.path.basename(x).split('_')[-1].replace('.mol2', '')))
        
        for mol2_file in mol2_files:
            pdbqt_file = os.path.join(folder, os.path.basename(mol2_file).replace('.mol2', '.pdbqt'))
            result = subprocess.run(['./lib/obabel', '-imol2', mol2_file, '-opdbqt', '-O', pdbqt_file],
                                    capture_output=True, text=True)
            
            self.contator += 1
            yield pdbqt_file

        for mol2_file in mol2_files:
            os.remove(mol2_file)
        
    def Maximum(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            self.maxim = len(lines)
