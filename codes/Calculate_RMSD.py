import subprocess
import os
from PySide6.QtCore import QThread, Signal


class Conversions(QThread):

    """
    A class responsible for performing RMSD calculations between only two 
    molecules using DockRMSD. It also creates two mol2 files in case the user requires them.
    """

    progress = Signal(int)
    finished = Signal()


    def __init__(self, file_1, file_2, folder_text):
        
        """
        Initializes the Conversions class with parameters for conversion.

        Parameters:
        -----------
        file_1 : str
            Path or name of the first file to be processed.
        file_2 : str
            Path or name of the second file to be processed.
        folder_text : str
            The destination folder path where the converted files will be saved.
        """

        super().__init__()
        self._running = True
        self.maxim = 0
        self.contator = 0
        self.file_1 = file_1
        self.file_2 = file_2
        self.des_folder = folder_text + "/"

    def run(self):

        # Create the destination folder
        os.makedirs(self.des_folder, exist_ok=True)

        # Create the names of the new files
        new_1 = self.des_folder + os.path.basename(self.file_1).split(".")[0] + ".mol2"
        new_2 = self.des_folder + os.path.basename(self.file_2).split(".")[0] + ".mol2"

        # Create the comand of the file 1 conversion to mol2 using openbabel
        command_1 = [
        'obabel', 
        '-ipdbqt', self.file_1, 
        '-omol2', "-O", new_1
        ]

        # Create the comand of the file 2 conversion to mol2 using openbabel
        command_2 = [
        'obabel', 
        '-ipdbqt', self.file_2, 
        '-omol2', "-O", new_2
        ]

        # Run the previous comands
        result1 = subprocess.run(command_1, capture_output=True, text=True)
        result2 = subprocess.run(command_2, capture_output=True, text=True)

        # Calcule the RMSD using DockRMSD for the two mol2 files
        result3 = subprocess.run(['./lib/DockRMSD', new_1, new_2], capture_output=True, text=True)

        # Write the results in a file called Results.txt
        with open(self.des_folder + "Results.txt", "w") as f:
            f.write(result3.stdout)
        
        self.contator +=1 
        self.progress.emit(self.contator)
        self.finished.emit()
    
    def Maximum(self):
        self.maxim = 1
    
    def stop(self):
        self._running = False
