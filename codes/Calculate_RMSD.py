import subprocess
import os
from PySide6.QtCore import QThread, Signal


class Conversions(QThread):
    progress = Signal(int)  # Señal para actualizar la barra de progreso
    finished = Signal()  # Señal para indicar que la conversión ha terminado


    def __init__(self, file_1, file_2, folder_text):
        super().__init__()
        self._running = True
        self.maxim = 0
        self.contator = 0
        self.file_1 = file_1
        self.file_2 = file_2
        self.des_folder = folder_text + "/"

    def run(self):
        os.makedirs(self.des_folder, exist_ok=True)

        new_1 = self.des_folder + os.path.basename(self.file_1).split(".")[0] + ".mol2"
        new_2 = self.des_folder + os.path.basename(self.file_2).split(".")[0] + ".mol2"

        command_1 = [
        './lib/obabel', 
        '-ipdbqt', self.file_1, 
        '-omol2', "-O", new_1
        ]

        command_2 = [
        './lib/obabel', 
        '-ipdbqt', self.file_2, 
        '-omol2', "-O", new_2
        ]

        result1 = subprocess.run(command_1, capture_output=True, text=True)
        result2 = subprocess.run(command_2, capture_output=True, text=True)

        result3 = subprocess.run(['./lib/DockRMSD', new_1, new_2], capture_output=True, text=True)

        with open(self.des_folder + "Results.txt", "w") as f:
            f.write(result3.stdout)
        
        self.contator +=1 

        self.progress.emit(self.contator)

        self.finished.emit()
    
    def Maximum(self):
        self.maxim = 1
    
    def stop(self):
        """Método para detener el hilo de forma segura"""
        self._running = False
