import os
import glob
import pandas as pd
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):
    progress = Signal(int)  # Señal para actualizar la barra de progreso
    finished = Signal()  # Señal para indicar que la conversión ha terminado

    def __init__(self, Number, logs, folder_text, index):
        super().__init__()
        self._running = True
        self.maxim = 0
        self.contator = 0
        self.Number = Number
        self.logs = logs + "/"
        self.des_folder = folder_text + "/"
        self.index = index

    
    def run(self):
        os.makedirs(self.des_folder, exist_ok=True)
        logs_files = glob.glob(self.logs + '*.log')

        final_list = []
        
        for log in logs_files:
            if not self._running:
                break
            temp_list = []

            with open(log,'r') as flog:
                lines = flog.readlines()
            
            for i, line in enumerate(lines):
                
                if i > 24 and i < len(lines):
                    line = os.path.basename(log).split(".")[0] + "       " + line
                    temp_list.append(line)
            
            final_list.append(temp_list)

            self.contator +=1

            self.progress.emit(self.contator)

        name = []
        pose = []
        score = []

        for main_line in final_list:
            if not self._running:
                break
            for line in main_line:
                name.append(line.split()[0])
                pose.append(line.split()[1])
                score.append(line.split()[2])

        my_dict = {}
        my_dict["name"] = name
        my_dict["pose"] = pose
        my_dict["score"] = score


        df = pd.DataFrame(my_dict)
        df['score'] = pd.to_numeric(df['score'], errors='coerce')
        df.sort_values(by="score", inplace=True, ascending=False)
        
        if self.index == 0:
            df.sort_values(by="score", inplace=True, ascending=True)
            df = df.head(int(self.Number)).copy()
        elif self.index == 1:
            df = df[df["score"] < float(self.Number)].copy()
            df.sort_values(by="score", inplace=True, ascending=True)

        df.to_csv(self.des_folder + "my_df.csv", index=False)

        self.finished.emit()
        

    def Maximum(self, logs):
        logs = logs + "/"
        pdbs_files = glob.glob(logs + '*.log')
        self.maxim = len(pdbs_files)
    
    def stop(self):
        """Método para detener el hilo de forma segura"""
        self._running = False
