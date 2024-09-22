import subprocess
import os
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):
    progress = Signal(int)  # Señal para actualizar la barra de progreso
    finished = Signal()  # Señal para indicar que la conversión ha terminado

    def __init__(self, file_path, folder_text):
        super().__init__()
        self.file_path = file_path
        self.folder_text = folder_text
        self.contator = 0
        self.maxim = 0
        self._is_running = True  # Variable para controlar la ejecución

    def run(self):
        real_folder = os.path.join(self.folder_text)
        os.makedirs(real_folder, exist_ok=True)

        # Crear listas vacías
        list_smiles = []
        list_ID = []

        # Guardar ID's y SMILES
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        for line in lines:
            parts = line.split()
            list_smiles.append(parts[0])
            list_ID.append(parts[1])

        self.maxim = len(list_smiles)  # Actualizar el máximo al iniciar

        for i, smile in enumerate(list_smiles):
            if not self._is_running:  # Verificar si se debe detener
                break
            try:
                subprocess.run(["./lib/obabel", "-:" + smile, "-omol2", "-O", os.path.join(real_folder, f'ligand_{i}.mol2'), '--gen3d'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                sdf = os.path.join(real_folder, f'ligand_{i}.mol2')
                pdbqt = os.path.join(real_folder, f'ligand_{i}.pdbqt')
                subprocess.run(['./lib/obabel', '-imol2', sdf, '-opdbqt', '-O', pdbqt], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                
                self.contator += 1
                os.remove(sdf)
                os.rename(pdbqt, os.path.join(real_folder, f"{list_ID[i]}.pdbqt"))
                
                # Emitir progreso
                self.progress.emit(self.contator)

            except Exception as e:
                print(f"Error processing {smile}: {str(e)}")

        self.finished.emit()  # Emitir señal al finalizar

    def stop(self):
        """Método para detener la ejecución del hilo."""
        self._is_running = False

    def Maximum(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            self.maxim = len(lines)