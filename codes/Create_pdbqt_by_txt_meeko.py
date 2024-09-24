import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):
    progress = Signal(int)

    """
    A class responsible for generating files in pdbqt format from a 
    single file in txt format containing several molecules in smiles format without their
    respective ID using meeko
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
        self.file_path = file_path
        self.folder_text = folder_text
        self.errors = {}
        self.maximum = 0
        self._running = True

    def run(self):
        list_smiles = []

        # Open the txt file
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        # save each smile on a list
        for line in lines:
            list_smiles.append(line.split()[0])

        self.maximum = len(list_smiles)

        # Make the destination folder
        os.makedirs(self.folder_text, exist_ok=True)

        for index, smile in enumerate(list_smiles):
            if not self._running:
                break

            name = str(index + 1)
            try:
                # Convert each smile on a mol rdkit object
                mol = Chem.MolFromSmiles(smile)
                if mol is None:
                    raise ValueError(f"Failed to create molecule from SMILES for index {name}")

                # Add hydrogens and optimize the molecule
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol)

                # Prepare the molecule
                preparator = MoleculePreparation()
                mol_setups = preparator.prepare(mol)
                if not mol_setups:
                    raise ValueError(f"No setups generated for index {name}")

                pdbqt_string = ""

                # Create the PDBQT for each molecule
                for setup in mol_setups:
                    result = PDBQTWriterLegacy.write_string(setup)
                    if isinstance(result, tuple):
                        pdbqt_string += result[0]
                    else:
                        pdbqt_string += result

                # Save each pbdqt on a file
                pdbqt_path = os.path.join(self.folder_text, name + ".pdbqt")
                with open(pdbqt_path, 'w') as archivo:
                    archivo.write(pdbqt_string)

            except Exception as e:
                self.errors[name] = str(e)

            self.progress.emit(index + 1)

        # Save the errors on a file
        if self.errors:
            error_file_path = os.path.join(self.folder_text, "errors.txt")
            with open(error_file_path, "w") as archivo:
                for key, value in self.errors.items():
                    archivo.write(f"{key}: {value}\n")

    def stop(self):
        self._running = False

    def Maximum(self, ruta: str):
        with open(ruta, 'r') as archivo:
            lineas = archivo.readlines()
        self.maximum = len(lineas)