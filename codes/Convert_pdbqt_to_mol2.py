import os
import glob

class Conversions():
    def __init__(self):
        self.maxim = 0
        self.contator = 0

    def conversions(self, pdbqt_files, folder_text):
        des_folder = folder_text + "/"
        ligands = pdbqt_files + "/"
        os.makedirs(des_folder, exist_ok=True)
        ligands_files = glob.glob(ligands + '*.pdbqt')
        
        
        for file in ligands_files:

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
                
                pdbqt_name = des_folder + os.path.basename(file).split(".")[0] + f"_{key}.pdbqt"
                mol2_name = des_folder + os.path.basename(file).split(".")[0] + f"_{key}.mol2"
                
                with open(pdbqt_name,'w') as f2:
                    for line in list:
                        f2.write(line)

            self.contator +=1
            yield
        
    
    def Maximum(self, ini_folder):
        ini_folder = ini_folder + "/"
        pdbs_files = glob.glob(ini_folder + '*.pdbqt')
        self.maxim = len(pdbs_files)

            
                