import subprocess
import os
import glob


class Conversions():
    def __init__(self):
        self.maxim = 0
        self.contator = 0
    
    def conversions(self, protein_file, ligands_folder, folder_text, x, y, z, sx, sy, sz, score):
        
        des_folder = folder_text + "/"
        ligands = ligands_folder + "/"
        os.makedirs(des_folder, exist_ok=True)
        ligands_files = glob.glob(ligands + '*.pdbqt')

        if score == 0:
            score = "vina"
        elif score == 1:
            score = "vinardo"
        elif score == 2:
            score = "dkoes_fast"

        with open(des_folder + "config.txt", "w") as file:
            file.write("----------Configuration employeed----------\n")
            file.write(f"Scoring: {score}\n")
            file.write(f"Center_X: {x}\n")
            file.write(f"Center_y: {y}\n")
            file.write(f"Center_z: {z}\n")
            file.write(f"size_x: {sx}\n")
            file.write(f"size_y: {sy}\n")
            file.write(f"size_z: {sz}")

        for ligand in ligands_files:
            basename = os.path.basename(ligand).split(".")[0]

            command = [
                './lib/smina',
                '-r', protein_file,
                '-l', ligand,
                '-o', des_folder + basename + "_docked.pdbqt",
                '--log', des_folder + basename + "_docked.log",
                '--center_x', x,
                '--center_y', y,
                '--center_z', z,
                '--size_x', sx,
                '--size_y', sy,
                '--size_z', sz,
                '--scoring', score
                ]
            
            # Ejecuta el comando
            subprocess.run(command, capture_output=True, text=True)

            self.contator +=1 
            yield
        
    
    def Maximum(self, ligands_folder):
        ligands_folder = ligands_folder + "/"
        pdbs_files = glob.glob(ligands_folder + '*.pdbqt')
        self.maxim = len(pdbs_files)

    