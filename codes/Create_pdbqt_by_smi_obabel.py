import subprocess
import os
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):

    """
    A class responsible for generating files in pdbqt format from a 
    single file in smi format containing several molecules in smiles format with their
    respective ID using openbabel
    """

    progress = Signal(int)
    finished = Signal()

    def __init__(self, file_path, folder_text):

        """
        Initializes the Conversions class with parameters for conversion.

        Parameters:
        -----------
        file_path : str
            Path to files in sdf format for the conversion to pdbqt
        folder_text : str
            The destination folder path where the converted files will be saved.
        """

        super().__init__()
        self.file_path = file_path
        self.folder_text = folder_text
        self.contator = 0
        self.maxim = 0
        self._is_running = True

    def run(self):
        real_folder = os.path.join(self.folder_text)
        os.makedirs(real_folder, exist_ok=True)

        # Make empty lists
        list_smiles = []
        list_ID = []

        # Open the smi file
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        # Save the smiles and the ID's on their respective list
        for line in lines:
            parts = line.split()
            list_smiles.append(parts[0])
            list_ID.append(parts[1])

        self.maxim = len(list_smiles)

        for i, smile in enumerate(list_smiles):
            if not self._is_running: 
                break
            try:
                # Run the openbabel command and convert to mol2
                subprocess.run(["obabel", "-:" + smile, "-omol2", "-O", os.path.join(real_folder, f'ligand_{i}.mol2'), '--gen3d'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                sdf = os.path.join(real_folder, f'ligand_{i}.mol2')
                pdbqt = os.path.join(real_folder, f'ligand_{i}.pdbqt')

                # Run the openbabel command and convert to pdbqt
                subprocess.run(['obabel', '-imol2', sdf, '-opdbqt', '-O', pdbqt], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                
                self.contator += 1
                os.remove(sdf)
                os.rename(pdbqt, os.path.join(real_folder, f"{list_ID[i]}.pdbqt"))
                self.progress.emit(self.contator)

            except Exception as e:
                print(f"Error processing {smile}: {str(e)}")

        self.finished.emit()

    def stop(self):
        self._is_running = False

    def Maximum(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            self.maxim = len(lines)