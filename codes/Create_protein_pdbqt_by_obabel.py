import subprocess
import os
import glob

class Conversions():
    def __init__(self):
        self.contator = 0
        self.maxim = 0
    
    def conversions(self, ini_folder, des_folder):
        ini_folder = ini_folder + "/"
        des_folder = des_folder + "/"
        os.makedirs(des_folder, exist_ok=True)
        pdbs_files = glob.glob(ini_folder + '*.pdb')

        for pdb in pdbs_files:
            final_name = des_folder + os.path.basename(os.path.splitext(pdb)[0]) + ".pdbqt"
            result = subprocess.run(["./lib/obabel", "-i", "pdb", pdb, "-xr", "-opdbqt", "-O", final_name], capture_output=True, text=True)
            self.contator += 1
            yield result
        
    def Maximum(self, ini_folder):
        ini_folder = ini_folder + "/"
        pdbs_files = glob.glob(ini_folder + '*.pdb')
        self.maxim = len(pdbs_files)



