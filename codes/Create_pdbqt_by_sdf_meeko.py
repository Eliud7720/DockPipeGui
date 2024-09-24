import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from PySide6.QtCore import QThread, Signal


class Conversions(QThread):

    """
    A class responsible for generating files in pdbqt format from a 
    single file in sdf format containing several molecules in that format
    using the Meeko library
    """

    progress = Signal(int)

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
        self.maxim = 0
        self._running = True

    def run(self):
        
        # Create the molecules list from the SDF file
        list_molecules = []
        supplier = Chem.SDMolSupplier(self.file_path)
        for mol in supplier:
            if mol is not None:
                list_molecules.append(mol)

        self.maxim = len(list_molecules)

        # Make the destination folder
        os.makedirs(self.folder_text, exist_ok=True)

        for index, mol in enumerate(list_molecules):
            if not self._running:
                break

            name = str(index + 1)

            try:
                # Add the hydrogens and optimize the molecule
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol)

                # Prepare the molecule and convert to PDBQT
                preparator = MoleculePreparation()
                mol_setups = preparator.prepare(mol)
                if not mol_setups:
                    raise ValueError(f"No setups generated for index {name}")

                pdbqt_string = ""
                for setup in mol_setups:
                    # Check if the return is a tuple
                    result = PDBQTWriterLegacy.write_string(setup)
                    if isinstance(result, tuple):
                        pdbqt_string += result[0]
                    else:
                        pdbqt_string += result

                # Save the pdbqt file
                pdbqt_path = os.path.join(self.folder_text, name + ".pdbqt")
                with open(pdbqt_path, 'w') as archivo:
                    archivo.write(pdbqt_string)

            except Exception as e:
                self.errors[name] = str(e)

            self.progress.emit(index + 1)

        # Save the errors on a file, if there are any
        if self.errors:
            error_file_path = os.path.join(self.folder_text, "errors.txt")
            with open(error_file_path, "w") as archivo:
                for key, value in self.errors.items():
                    archivo.write(f"{key}: {value}\n")

    def stop(self):
        self._running = False

    def Maximum(self, ruta: str):
        list_molecules = []
        supplier = Chem.SDMolSupplier(ruta)
        for mol in supplier:
            if mol is not None:
                list_molecules.append(mol)

        self.maxim = len(list_molecules)
