import subprocess
import os
from rdkit import Chem
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):
    progress = Signal(int)  # Señal para actualizar la barra de progreso
    finished = Signal()  # Señal para indicar que la conversión ha terminado

    def __init__(self, file, folder):
        super().__init__()
        self.contator = 0
        self.maxim = 0
        self.folder = folder
        self.file = file
        self._is_running = True  # Variable para controlar la ejecución

    def run(self):
        real_folder = os.path.join(self.folder)
        os.makedirs(real_folder, exist_ok=True)

        # Lista de SMILES
        smiles_list = []
        supplier = Chem.SDMolSupplier(self.file)

        for mol in supplier:
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                smiles_list.append(smiles)

        self.maxim = len(smiles_list)  # Establecer el máximo

        for i, smile in enumerate(smiles_list):
            if not self._is_running:  # Verificar si se debe detener
                break
            
            try:
                subprocess.run(
                    ["./lib/obabel", "-:" + smile, "-omol2", "-O", os.path.join(real_folder, f'ligand_{i}.mol2'), '--gen3d'],
                    check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
                )
                sdf = os.path.join(real_folder, f'ligand_{i}.mol2')
                pdbqt = os.path.join(real_folder, f'ligand_{i}.pdbqt')
                subprocess.run(
                    ['./lib/obabel', '-imol2', sdf, '-opdbqt', '-O', pdbqt],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )
                self.contator += 1
                os.remove(sdf)

                # Emitir progreso
                self.progress.emit(self.contator)

            except Exception as e:
                print(f"Error processing {smile}: {str(e)}")

        self.finished.emit()  # Emitir señal al finalizar

    def stop(self):
        """Método para detener la ejecución del hilo."""
        self._is_running = False

    def Maximum(self, file):
        list_molecules = []

        supplier = Chem.SDMolSupplier(file)
        for mol in supplier:
            if mol is not None:
                list_molecules.append(mol)

        self.maxim = len(list_molecules)
