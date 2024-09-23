import subprocess
import os
import numpy as np
from PySide6.QtCore import QThread, Signal


class Conversions(QThread):
    progress = Signal(int)  # Señal para actualizar la barra de progreso
    finished = Signal()  # Señal para indicar que la conversión ha terminado

    def __init__(self, protein_file, ligands_folder, folder_text, score):
        self.maxim = 0
        self.contator = 0
        super().__init__()
        self.des_folder = folder_text + "/"
        self.ligand = ligands_folder
        self.protein_file = protein_file
        self.score = score
        self._running = True
    
    def run(self):
        
        os.makedirs(self.des_folder, exist_ok=True)

        
        if self.score == 0:
            self.score = "vina"
        elif self.score == 1:
            self.score = "vinardo"
        elif self.score == 2:
            self.score = "dkoes"
        
        xyz_list = []
        
        
        with open(self.ligand, 'r') as file:
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

            if not self._running:  # Verificar si se debe detener el hilo
                break

            element = np.array(element, dtype = float)

            new_distance = np.linalg.norm(element-geometric_array)

            if new_distance > distance:
                distance = new_distance
        
        size_X = distance + 7
        size_Y = distance + 7
        size_Z = distance + 7

        

        basename = os.path.basename(self.ligand).split(".")[0]

        with open(self.des_folder + "config.txt", "w") as file:
            file.write("----------Configuration employeed----------\n")
            file.write(f"Scoring: {self.score}\n")
            file.write(f"Center_X: {round(geometric_center[0], 2)}\n")
            file.write(f"Center_y: {round(geometric_center[1], 2)}\n")
            file.write(f"Center_z: {round(geometric_center[2], 2)}\n")
            file.write(f"size_x: {size_X}\n")
            file.write(f"size_y: {size_Y}\n")
            file.write(f"size_z: {size_Z}")

        command = [
            './lib/smina',
            '-r', self.protein_file,
            '-l', self.ligand,
            '-o', self.des_folder + basename + "_docked.pdbqt",
            '--log', self.des_folder + basename + "_docked.log",
            '--center_x', str(round(geometric_center[0], 2)),
            '--center_y', str(round(geometric_center[1], 2)),
            '--center_z', str(round(geometric_center[2], 2)),
            '--size_x', str(size_X),
            '--size_y', str(size_Y),
            '--size_z', str(size_Z),
            '--scoring', self.score
            ]
        
        subprocess.run(command, capture_output=True, text=True)

        self.progress.emit(self.contator + 1)

        self.finished.emit()

    def Maximum(self, ligands_folder):
        self.maxim = 1
    
    def stop(self):
        """Método para detener el hilo de forma segura"""
        self._running = False

    