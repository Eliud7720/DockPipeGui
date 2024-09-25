import subprocess
import os
from rdkit import Chem
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):

    """
    A class responsible for generating files in pdbqt format from a 
    single file in sdf format containing several molecules in that format
    using the openbabel software
    """

    progress = Signal(int)
    finished = Signal()

    def __init__(self, file, folder):

        """
        Initializes the Conversions class with parameters for conversion.

        Parameters:
        -----------
        file : str
            Path to files in sdf format for the conversion to pdbqt
        folder : str
            The destination folder path where the converted files will be saved.
        """

        super().__init__()
        self.contator = 0
        self.maxim = 0
        self.folder = folder
        self.file = file
        self._is_running = True

    def run(self):

        # Make the destination folder
        real_folder = os.path.join(self.folder)
        os.makedirs(real_folder, exist_ok=True)

        # Smile list
        smiles_list = []
        supplier = Chem.SDMolSupplier(self.file)

        # Convert to mol object each smiles
        for mol in supplier:
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                smiles_list.append(smiles)

        self.maxim = len(smiles_list)  # Establish the maximum

        for i, smile in enumerate(smiles_list):
            if not self._is_running:
                break
            
            try:
                # Run the conversion to mol2 to each molecule
                subprocess.run(
                    ["obabel", "-:" + smile, "-omol2", "-O", os.path.join(real_folder, f'ligand_{i}.mol2'), '--gen3d'],
                    check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
                )
                sdf = os.path.join(real_folder, f'ligand_{i}.mol2')
                pdbqt = os.path.join(real_folder, f'ligand_{i}.pdbqt')

                # Run the conversion to pdbqt to each molecule
                subprocess.run(
                    ['obabel', '-imol2', sdf, '-opdbqt', '-O', pdbqt],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )
                self.contator += 1
                os.remove(sdf)

                self.progress.emit(self.contator)

            except Exception as e:
                print(f"Error processing {smile}: {str(e)}")

        self.finished.emit()

    def stop(self):
        self._is_running = False

    def Maximum(self, file):
        list_molecules = []

        supplier = Chem.SDMolSupplier(file)
        for mol in supplier:
            if mol is not None:
                list_molecules.append(mol)

        self.maxim = len(list_molecules)
