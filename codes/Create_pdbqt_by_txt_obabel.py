import subprocess
import os
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):
    progress = Signal(int)
    finished = Signal()

    """
    A class responsible for generating files in pdbqt format from a 
    single file in txt format containing several molecules in smiles format without their
    respective ID using openbabel
    """

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
        self.contator = 0
        self.maxim = 0
        self._running = True
        self.file_path = file_path
        self.folder_text = folder_text

    def run(self):
        real_folder = self.folder_text + "/"
        os.makedirs(real_folder, exist_ok=True)

        # Create the smiles list
        list_smiles = []

        # open the smiles file
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        # Save the smiles on a list
        for line in lines:
            list_smiles.append(line.strip())

        self.maxim = len(list_smiles)

        for i, smile in enumerate(list_smiles):
            if not self._running:
                break
            
            try:
                # Run the command to convert each smiles to mol2
                subprocess.run(
                    ["obabel", "-:" + smile, "-omol2", "-O", os.path.join(real_folder, f'ligand_{i}.mol2'), '--gen3d'],
                    check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
                )
                sdf = os.path.join(real_folder, f'ligand_{i}.mol2')
                pdbqt = os.path.join(real_folder, f'ligand_{i}.pdbqt')

                # Run the command to convert each smiles to pdbqt
                subprocess.run(
                    ['obabel', '-imol2', sdf, '-opdbqt', '-O', pdbqt],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )
                self.contator += 1
                os.remove(sdf)

                self.progress.emit(self.contator)

            except subprocess.CalledProcessError as e:
                print(f"Error processing {smile}: {str(e)}")

        self.finished.emit()

    def stop(self):
        self._running = False

    def Maximum(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            self.maxim = len(lines)
