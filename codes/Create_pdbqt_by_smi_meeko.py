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
    single file in smi format containing several molecules in smiles format with their
    respective ID using the Meeko library
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
        self.contator = 0
        self.errors = {}
        self.maxim = 0
        self._running = True

    def run(self):
        """
        Main function responsible for performing conversions.
        This function is executed when the thread starts.
        """

        # Make of the lists
        list_smiles = []
        list_ID = []

        # Open the smi file
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        # Save each smiles and the ID on their respective list
        for line in lines:
            list_smiles.append(line.split()[0])
            list_ID.append(line.split()[1])

        self.maxim = len(list_smiles)

        # Make the destination directory
        os.makedirs(self.folder_text, exist_ok=True)

        for index, smile in enumerate(list_smiles):
            if not self._running:
                break

            name = list_ID[index]
            try:
                # Convert each smiles to mol rdkit object 
                mol = Chem.MolFromSmiles(smile)
                if mol is None:
                    raise ValueError(f"Failed to create molecule from SMILES for ID {name}")

                # Add hydrogends and optimize the molecule
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol)

                # Prepare the molecule object
                preparator = MoleculePreparation()
                mol_setups = preparator.prepare(mol)
                if not mol_setups:
                    raise ValueError(f"No setups generated for ID {name}")

                # Calculate the pdbqt for each molecule
                pdbqt_string = ""
                for setup in mol_setups:
                    result = PDBQTWriterLegacy.write_string(setup)
                    if isinstance(result, tuple):
                        pdbqt_string += result[0]
                    else:
                        pdbqt_string += result

                # Write the pdbqt files
                pdbqt_path = os.path.join(self.folder_text, name + ".pdbqt")
                with open(pdbqt_path, 'w') as file:
                    file.write(pdbqt_string)

            except Exception as e:
                self.errors[name] = str(e)

            self.progress.emit(index + 1)

        if self.errors:
            error_file_path = os.path.join(self.folder_text, "errors.txt")
            with open(error_file_path, "w") as file:
                for key, value in self.errors.items():
                    file.write(f"{key}: {value}\n")
    
    def stop(self):
        self._running = False
    
    def Maximum(self, path: str):
        list_smiles = []
        list_ID = []

        with open(path, 'r') as file:
            lineas = file.readlines()

        for linea in lineas:
            list_smiles.append(linea.split()[0])
            list_ID.append(linea.split()[1])

        self.maxim = len(list_smiles)
