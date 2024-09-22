import subprocess
import os
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):  # Hereda de QThread
    progress = Signal(int)  # Señal para actualizar la barra de progreso
    finished = Signal()  # Señal para indicar que la conversión ha terminado

    def __init__(self, file_path, folder_text):
        super().__init__()  # Inicializar la clase base
        self.contator = 0
        self.maxim = 0
        self._running = True
        self.file_path = file_path
        self.folder_text = folder_text

    def run(self):
        real_folder = self.folder_text + "/"
        os.makedirs(real_folder, exist_ok=True)

        # Crear lista de SMILES
        list_smiles = []

        # Guardar SMILES
        with open(self.file_path, 'r') as file:
            lines = file.readlines()

        for line in lines:
            list_smiles.append(line.strip())  # Asegúrate de eliminar saltos de línea

        self.maxim = len(list_smiles)  # Establecer el máximo

        for i, smile in enumerate(list_smiles):
            if not self._running:
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
                os.remove(sdf)  # Eliminar el archivo .mol2

                # Emitir progreso
                self.progress.emit(self.contator)

            except subprocess.CalledProcessError as e:
                print(f"Error processing {smile}: {str(e)}")

        self.finished.emit()  # Emitir señal al finalizar

    def stop(self):
        """Método para detener la ejecución del hilo."""
        self._running = False

    def Maximum(self, file):
        with open(file, 'r') as f:
            lines = f.readlines()
            self.maxim = len(lines)
