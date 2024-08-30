import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy


class Conversiones():
    """
    Función encargada de realizar conversiones a pdbqt y sdf a partir 
    de un documento de texto lleno de SMILES con su respectivo identificador
    """

    def __init__(self):
        self.contador = 0
        self.errores = {}
        self.maximo = 0

    def conversion(self, ruta: str):
        """
        Funcion principal de la función encargada de realizar las conversiones.

        ruta: Seleccionar la ruta del archivo.
        """
        # Creación de las listas
        list_smiles = []
        list_ID = []

        # Abrir archivo
        with open(ruta, 'r') as archivo:
            lineas = archivo.readlines()

        # Guardar en cada lista los SMILES y los ID
        for linea in lineas:
            list_smiles.append(linea.split()[0])
            list_ID.append(linea.split()[1])

        self.maximo = len(list_smiles)

        # Crear una lista de objetos mol con la lista de smiles
        list_mol = []

        for smile in list_smiles:
            list_mol.append(Chem.MolFromSmiles(smile))

        # Progreso de la conversión
        progress = []

        # Crear la carpeta "PDBQT files" si no existe
        os.makedirs("PDBQT files", exist_ok=True)

        # Crear un archivo sdf para cada molécula
        for molec in list_mol:
            try:
                mol = Chem.AddHs(molec)
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol)

                # Guardar la molécula en un archivo SDF en la carpeta "PDBQT files"
                nombre = list_ID[self.contador]
                sdf_path = os.path.join("PDBQT files", nombre + ".sdf")
                with Chem.SDWriter(sdf_path) as w:
                    w.write(mol)

            except Exception as e:
                self.errores[nombre] = e

            self.contador += 1

            # Actualizar el progreso de la conversión
            progress.append(self.contador)

            # Devolver el progreso actual
            yield self.contador, self.maximo, progress

        # Convertir SDF a PDBQT y guardar en la carpeta "PDBQT files"
        for id in list_ID:
            try:
                sdf_path = os.path.join("PDBQT files", id + ".sdf")

                # there is one molecule in this SD file, this loop iterates just once
                for mol in Chem.SDMolSupplier(sdf_path, removeHs=False):
                    preparator = MoleculePreparation()
                    mol_setups = preparator.prepare(mol)
                    for setup in mol_setups:
                        pdbqt_string = PDBQTWriterLegacy.write_string(setup)

                pdbqt_string = pdbqt_string[0]

                pdbqt_path = os.path.join("PDBQT files", id + ".pdbqt")
                with open(pdbqt_path, 'a') as archivo:
                    # Escribe más contenido en el archivo
                    archivo.write(pdbqt_string)

            except Exception as e:
                self.errores[id] = e

        # Guardar los errores en un archivo si hay alguno
        if self.errores:
            with open("errores.txt", "w") as archivo:
                for key, value in self.errores.items():
                    archivo.write(f"{key}: {value}\n")


class Maximo():
    def __init__(self):
        self.contador = 0
        self.errores = {}
        self.maximo = 0

    def contar_maximo(self, ruta: str):
        # Creación de las listas
        list_smiles = []
        list_ID = []

        # Abrir archivo
        with open(ruta, 'r') as archivo:
            lineas = archivo.readlines()

        # Guardar en cada lista los SMILES y los ID
        for linea in lineas:
            list_smiles.append(linea.split()[0])
            list_ID.append(linea.split()[1])

        self.maximo = len(list_smiles)
