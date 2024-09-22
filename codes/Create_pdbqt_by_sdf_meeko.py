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
        self.maxim = 0
        self._running = True  # Variable para controlar el hilo

    def run(self):
        """
        Función principal responsable de realizar las conversiones.
        """
        # Crear lista de moléculas desde el archivo SDF
        list_molecules = []
        supplier = Chem.SDMolSupplier(self.file_path)
        for mol in supplier:
            if mol is not None:
                list_molecules.append(mol)

        self.maxim = len(list_molecules)

        # Crear la carpeta si no existe
        os.makedirs(self.folder_text, exist_ok=True)

        for index, mol in enumerate(list_molecules):
            if not self._running:  # Verificar si se debe detener el hilo
                break

            name = str(index + 1)  # Usar el índice como nombre de archivo

            try:
                # Agregar hidrógenos, embeber y optimizar la molécula
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
                    # Verificar si el retorno es una tupla
                    result = PDBQTWriterLegacy.write_string(setup)
                    if isinstance(result, tuple):
                        pdbqt_string += result[0]
                    else:
                        pdbqt_string += result

                # Guardar el archivo PDBQT
                pdbqt_path = os.path.join(self.folder_text, name + ".pdbqt")
                with open(pdbqt_path, 'w') as archivo:
                    archivo.write(pdbqt_string)

            except Exception as e:
                self.errors[name] = str(e)

            # Emitir progreso
            self.progress.emit(index + 1)

        # Guardar errores en un archivo, si los hay
        if self.errors:
            error_file_path = os.path.join(self.folder_text, "errors.txt")
            with open(error_file_path, "w") as archivo:
                for key, value in self.errors.items():
                    archivo.write(f"{key}: {value}\n")

    def stop(self):
        """Método para detener el hilo de forma segura"""
        self._running = False

    def Maximum(self, ruta: str):
        """Función para contar el número máximo de moléculas en el archivo SDF"""
        list_molecules = []
        supplier = Chem.SDMolSupplier(ruta)
        for mol in supplier:
            if mol is not None:
                list_molecules.append(mol)

        self.maxim = len(list_molecules)
