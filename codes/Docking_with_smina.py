import subprocess
import os
import glob
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):

    """
    A class responsible for carrying out molecular docking using smina in a simpler way
    """

    progress = Signal(int)
    finished = Signal()

    def __init__(self, protein_file, ligands_folder, folder_text, x, y, z, sx, sy, sz, score):
        
        """
        Initializes the Conversions class with parameters for conversion.

        Parameters:
        -----------
        protein_file : str
            Path to the file of the protein in pdbqt format
        ligands_folder : str
            Path to the ligands folder in pdbqt format
        folder_text: str
            The destination folder path where the converted files will be saved.
        x : float
            X center coordinates for docking
        y : float
            Y center coordinates for docking
        z : float
            Z center coordinates for docking
        sx : float
            x size of the box
        sy : float
            y size of the box
        sz : float
            z size of the box
        score : int
            The method of scoring (vina, vinardo or dkoes)
        """

        super().__init__()
        self.maxim = 0
        self.contator = 0
        self.protein_file = protein_file
        self.ligands = ligands_folder + "/"
        self.des_folder = folder_text + "/"
        self.x = x
        self.y = y
        self.z = z
        self.sx = sx
        self.sy = sy
        self.sz = sz
        self.score = score
        self._running = True
    
    def run(self):
        
        # Make the destination folder
        os.makedirs(self.des_folder, exist_ok=True)
        ligands_files = glob.glob(self.ligands + '*.pdbqt')

        if self.score == 0:
            self.score = "vina"
        elif self.score == 1:
            self.score = "vinardo"
        elif self.score == 2:
            self.score = "dkoes_fast"

        # Create a file with the data used in case the user needs it
        with open(self.des_folder + "config.txt", "w") as file:
            file.write("----------Configuration employeed----------\n")
            file.write(f"Scoring: {self.score}\n")
            file.write(f"Center_X: {self.x}\n")
            file.write(f"Center_y: {self.y}\n")
            file.write(f"Center_z: {self.z}\n")
            file.write(f"size_x: {self.sx}\n")
            file.write(f"size_y: {self.sy}\n")
            file.write(f"size_z: {self.sz}")

        for ligand in ligands_files:
            if not self._running:
                break
            
            # Recover the basename of the ligand
            basename = os.path.basename(ligand).split(".")[0]

            # Command to run
            command = [
                'smina',
                '-r', self.protein_file,
                '-l', ligand,
                '-o', self.des_folder + basename + "_docked.pdbqt",
                '--log', self.des_folder + basename + "_docked.log",
                '--center_x', self.x,
                '--center_y', self.y,
                '--center_z', self.z,
                '--size_x', self.sx,
                '--size_y', self.sy,
                '--size_z', self.sz,
                '--scoring', self.score
                ]
            
            # Run the command
            subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            self.contator +=1 
            self.progress.emit(self.contator)
        
        self.finished.emit()
        
    def Maximum(self, ligands_folder):
        ligands_folder = ligands_folder + "/"
        pdbs_files = glob.glob(ligands_folder + '*.pdbqt')
        self.maxim = len(pdbs_files)
    
    def stop(self):
        self._running = False

    