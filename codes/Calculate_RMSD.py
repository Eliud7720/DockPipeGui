import subprocess
import os
import glob

class Conversions():
    def __init__(self):
        self.maxim = 0
        self.contator = 0

    def conversions(self, file_1, file_2, folder_text):
        des_folder = folder_text + "/"
        os.makedirs(des_folder, exist_ok=True)

        new_1 = des_folder + os.path.basename(file_1).split(".")[0] + ".sdf"
        new_2 = des_folder + os.path.basename(file_2).split(".")[0] + ".sdf"

        command_1 = [
        './lib/obabel', 
        '-ipdbqt', file_1, 
        '-omol2', "-O", new_1
        ]

        command_2 = [
        './lib/obabel', 
        '-ipdbqt', file_2, 
        '-omol2', "-O", new_2
        ]

        result1 = subprocess.run(command_1, capture_output=True, text=True)
        result2 = subprocess.run(command_2, capture_output=True, text=True)

        result3 = subprocess.run(['./lib/DockRMSD', new_1, new_2], capture_output=True, text=True)

        with open(des_folder + "Results.txt", "w") as f:
            f.write(result3.stdout)


    
    def Maximum(self,):
        self.maxim = 1