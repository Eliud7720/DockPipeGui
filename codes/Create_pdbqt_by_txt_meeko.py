import os
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):
    progress = Signal(int)  # Señal para actualizar la barra de progreso

    def __init__(self, file_path, folder_text):
        super().__init__()
        self.file_path = file_path
        self.folder_text = folder_text
        self.errors = {}
        self.maximum = 0
        self._running = True

    def run(self):
        # Crear la lista de SMILES
        list_smiles = []

        # Abrir el archivo
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        # Guardar SMILES en la lista
        for line in lines:
            list_smiles.append(line.split()[0])

        self.maximum = len(list_smiles)

        # Crear la carpeta si no existe
        os.makedirs(self.folder_text, exist_ok=True)

        for index, smile in enumerate(list_smiles):
            if not self._running:
                break  # Detener el hilo si se ha solicitado

            name = str(index + 1)
            try:
                # Convertir SMILES a molécula
                mol = Chem.MolFromSmiles(smile)
                if mol is None:
                    raise ValueError(f"Failed to create molecule from SMILES for index {name}")

                # Agregar hidrógenos, incrustar y optimizar la molécula
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.MMFFOptimizeMolecule(mol)

                # Preparar la molécula y convertir a PDBQT
                preparator = MoleculePreparation()
                mol_setups = preparator.prepare(mol)
                if not mol_setups:
                    raise ValueError(f"No setups generated for index {name}")

                pdbqt_string = ""
                for setup in mol_setups:
                    result = PDBQTWriterLegacy.write_string(setup)
                    if isinstance(result, tuple):
                        pdbqt_string += result[0]
                    else:
                        pdbqt_string += result

                # Guardar PDBQT en archivo
                pdbqt_path = os.path.join(self.folder_text, name + ".pdbqt")
                with open(pdbqt_path, 'w') as archivo:
                    archivo.write(pdbqt_string)

            except Exception as e:
                self.errors[name] = str(e)

            # Emitir progreso
            self.progress.emit(index + 1)

        # Guardar errores en un archivo
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