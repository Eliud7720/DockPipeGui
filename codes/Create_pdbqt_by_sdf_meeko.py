import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy

class Conversions:
    """
    Function responsible for performing conversions 
    to pdbqt from an SDF file containing multiple molecules.
    """

    def __init__(self):
        self.errors = {}
        self.maximum = 0

    def conversion(self, path: str, folder: str):
        """
        Main function responsible for performing conversions.

        path: Select the file path.
        """
        # Create the list of molecules
        list_molecules = []

        # Open the file and read molecules
        supplier = Chem.SDMolSupplier(path)
        for mol in supplier:
            if mol is not None:
                list_molecules.append(mol)
        
        self.maximum = len(list_molecules)

        # Create the directory if it does not exist
        os.makedirs(folder, exist_ok=True)

        for index, mol in enumerate(list_molecules):
            name = str(index + 1)  # Use the index as the file name (1, 2, 3, etc.)
            try:
                # Add hydrogens, embed, and optimize the molecule
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
                    # Check if the return type is a tuple
                    result = PDBQTWriterLegacy.write_string(setup)
                    if isinstance(result, tuple):
                        pdbqt_string += result[0]  # Assuming the first element of the tuple is the string
                    else:
                        pdbqt_string += result

                # Save PDBQT to file
                pdbqt_path = os.path.join(folder, name + ".pdbqt")
                with open(pdbqt_path, 'w') as archivo:
                    archivo.write(pdbqt_string)

            except Exception as e:
                self.errors[name] = str(e)

            # Yield progress
            yield index + 1, self.maximum, list_molecules[:index + 1]

        # Save errors in a file
        if self.errors:
            error_file_path = os.path.join(folder, "errors.txt")
            with open(error_file_path, "w") as archivo:
                for key, value in self.errors.items():
                    archivo.write(f"{key}: {value}\n")


class Maximum:
    def __init__(self):
        self.contador = 0
        self.errores = {}
        self.maximo = 0

    def contar_maximo(self, ruta: str):
        list_molecules = []

        supplier = Chem.SDMolSupplier(ruta)
        for mol in supplier:
            if mol is not None:
                list_molecules.append(mol)

        self.maximo = len(list_molecules)
