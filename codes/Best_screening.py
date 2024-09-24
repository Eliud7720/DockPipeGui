import os
import glob
import pandas as pd
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):

    """
    A class in charge of carrying out a screening and extraction of the 
    best molecules according to their docking score from files generated 
    by the molecular docking '.log'
    """

    progress = Signal(int)
    finished = Signal()

    def __init__(self, Number, logs, folder_text, index):

        """
        Initializes the Conversions class with parameters for conversion.

        Parameters:
        -----------
        Number : float
            Number of molecules requested by the user or requested interval
        logs : str
            The path to the logs directory where logs will be stored.
        folder_text : str
            The destination folder path where the converted files will be saved.
        index : int
            Represents the index of the combo box to know what type of operation the user will do.
        """
    
        super().__init__()
        self._running = True
        self.maxim = 0
        self.contator = 0
        self.Number = Number
        self.logs = logs + "/"
        self.des_folder = folder_text + "/"
        self.index = index

    
    def run(self):

        # Make the destination directory
        os.makedirs(self.des_folder, exist_ok=True)

        # Search the logs files 
        logs_files = glob.glob(self.logs + '*.log')
        final_list = []
        
        # Extract log lines from each log file
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

        # Create the name, pose and score lists
        name = []
        pose = []
        score = []

        # fill the lists with their respective information
        for main_line in final_list:
            if not self._running:
                break
            for line in main_line:
                name.append(line.split()[0])
                pose.append(line.split()[1])
                score.append(line.split()[2])

        # Create a dictionary with this information
        my_dict = {}
        my_dict["name"] = name
        my_dict["pose"] = pose
        my_dict["score"] = score

        # Create a dataframe with this information
        df = pd.DataFrame(my_dict)
        df['score'] = pd.to_numeric(df['score'], errors='coerce')
        df.sort_values(by="score", inplace=True, ascending=False)
        
        # Determine the best ligands according to the type required by the user
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
        self._running = False
