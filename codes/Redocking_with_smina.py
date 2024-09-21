import subprocess
import os
import numpy as np


class Conversions():
    def __init__(self):
        self.maxim = 0
        self.contator = 0
    
    def conversions(self, protein_file, ligands_folder, folder_text, score):
        
        des_folder = folder_text + "/"
        ligand = ligands_folder
        os.makedirs(des_folder, exist_ok=True)

        
        if score == 0:
            score = "vina"
        elif score == 1:
            score = "vinardo"
        elif score == 2:
            score = "dkoes"
        
        xyz_list = []
        
        
        with open(ligand, 'r') as file:
            lines = file.readlines()
        
        for i, line in enumerate(lines):
            if i > 3 and i < len(lines)-3:
                line = line[30:]
                xyz_list.append([line.split()[0], line.split()[1], line.split()[2]])
        
        array = np.array(xyz_list)
        data_array = np.array(array, dtype=float)

        geometric_center = np.mean(data_array, axis=0)
        geometric_array = np.array(geometric_center)

        distance = 0

        for element in xyz_list:
            element = np.array(element, dtype = float)

            new_distance = np.linalg.norm(element-geometric_array)

            if new_distance > distance:
                distance = new_distance
        
        size_X = distance + 7
        size_Y = distance + 7
        size_Z = distance + 7


        basename = os.path.basename(ligand).split(".")[0]

        with open(des_folder + "config.txt", "w") as file:
            file.write("----------Configuration employeed----------\n")
            file.write(f"Scoring: {score}\n")
            file.write(f"Center_X: {round(geometric_center[0], 2)}\n")
            file.write(f"Center_y: {round(geometric_center[1], 2)}\n")
            file.write(f"Center_z: {round(geometric_center[2], 2)}\n")
            file.write(f"size_x: {size_X}\n")
            file.write(f"size_y: {size_Y}\n")
            file.write(f"size_z: {size_Z}")

        command = [
            './lib/smina',
            '-r', protein_file,
            '-l', ligand,
            '-o', des_folder + basename + "_docked.pdbqt",
            '--log', des_folder + basename + "_docked.log",
            '--center_x', str(round(geometric_center[0], 2)),
            '--center_y', str(round(geometric_center[1], 2)),
            '--center_z', str(round(geometric_center[2], 2)),
            '--size_x', str(size_X),
            '--size_y', str(size_Y),
            '--size_z', str(size_Z),
            '--scoring', score
            ]
        
        subprocess.run(command, capture_output=True, text=True)

        self.contator +=1 


        
    
    def Maximum(self, ligands_folder):
        self.maxim = 1

    