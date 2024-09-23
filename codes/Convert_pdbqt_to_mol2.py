import os
import glob
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):
    progress = Signal(int)  # Señal para actualizar la barra de progreso
    finished = Signal()  # Señal para indicar que la conversión ha terminado

    def __init__(self, pdbqt_files, folder_text):
        super().__init__()
        self.maxim = 0
        self.contator = 0
        self.des_folder = folder_text + "/"
        self.ligands = pdbqt_files + "/"
        self._running = True

    def run(self):
        os.makedirs(self.des_folder, exist_ok=True)
        ligands_files = glob.glob(self.ligands + '*.pdbqt')
        
        
        for file in ligands_files:

            if not self._running:
                break

            with open(file, 'r') as f:
                lines = f.readlines()

            my_dict = {}
            temp_list = []
            con = 0

            for line in lines:

                if line.startswith("MODEL"):
                    con +=1
                    my_dict[con] = []
                
                my_dict[con].append(line)


            for key, list in my_dict.items():
                
                pdbqt_name = self.des_folder + os.path.basename(file).split(".")[0] + f"_{key}.pdbqt"
                mol2_name = self.des_folder + os.path.basename(file).split(".")[0] + f"_{key}.mol2"
                
                with open(pdbqt_name,'w') as f2:
                    for line in list:
                        f2.write(line)

            self.contator +=1
            self.progress.emit(self.contator)
        
        self.finished.emit()
        
    def Maximum(self, ini_folder):
        ini_folder = ini_folder + "/"
        pdbs_files = glob.glob(ini_folder + '*.pdbqt')
        self.maxim = len(pdbs_files)
    
    def stop(self):
        """Método para detener el hilo de forma segura"""
        self._running = False

            
                