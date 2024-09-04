import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy


class Conversions:
    """
    Function responsible for performing conversions 
    to pdbqt from a text document filled with SMILES with their respective identifier.
    """

    def __init__(self):
        self.contator = 0
        self.errors = {}
        self.maximum = 0

    def conversion(self, path: str, folder: str):
        """
        Main function responsible for performing conversions.

        path: Select the file path.
        str: The name of the folder.
        """
        # Creation of the lists
        list_smiles = []
        list_ID = []

        # Open the file
        with open(path, 'r') as file:
            lines = file.readlines()

        # Save SMILES and IDs in each list
        for line in lines:
            list_smiles.append(line.split()[0])
            list_ID.append(line.split()[1])

        self.maximum = len(list_smiles)

        # Create the folder if not exists
        os.makedirs(folder, exist_ok=True)

        for index, smile in enumerate(list_smiles):
            name = list_ID[index]
            try:
                # Convert SMILES to molecule
                mol = Chem.MolFromSmiles(smile)
                if mol is None:
                    raise ValueError(f"Failed to create molecule from SMILES for ID {name}")

                # Add hydrogens, embed, and optimize the molecule
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol)

                # Prepare the molecule and convert to PDBQT
                preparator = MoleculePreparation()
                mol_setups = preparator.prepare(mol)
                if not mol_setups:
                    raise ValueError(f"No setups generated for ID {name}")

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
            yield index + 1, self.maximum, list_ID[:index + 1]

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
        list_smiles = []
        list_ID = []

        with open(ruta, 'r') as archivo:
            lineas = archivo.readlines()

        for linea in lineas:
            list_smiles.append(linea.split()[0])
            list_ID.append(linea.split()[1])

        self.maximo = len(list_smiles)
