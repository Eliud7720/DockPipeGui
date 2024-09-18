import os
import glob
import pandas as pd

class Conversions():
    def __init__(self):
        self.maxim = 0
        self.contator = 0
    
    def conversions(self, Number, logs, folder_text, index):
        des_folder = folder_text + "/"
        logs = logs + "/"
        os.makedirs(des_folder, exist_ok=True)
        logs_files = glob.glob(logs + '*.log')

        final_list = []
        
        for log in logs_files:
            temp_list = []

            with open(log,'r') as flog:
                lines = flog.readlines()
            
            for i, line in enumerate(lines):
                
                if i > 24 and i < len(lines):
                    line = os.path.basename(log).split(".")[0] + "       " + line
                    temp_list.append(line)
            
            final_list.append(temp_list)

            self.contator +=1
        yield

        name = []
        pose = []
        score = []

        for main_line in final_list:
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
        
        if index == 0:
            df = df.head(int(Number)).copy()
        elif index == 1:
            df = df[df["score"] < float(Number)].copy()

        df.sort_values(by="score", inplace=True, ascending=True)
        df.to_csv(des_folder + "my_df.csv", index=False)
        





    def Maximum(self, logs):
        logs = logs + "/"
        pdbs_files = glob.glob(logs + '*.log')
        self.maxim = len(pdbs_files)
