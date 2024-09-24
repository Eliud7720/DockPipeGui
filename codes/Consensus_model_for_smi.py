import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from joblib import load
import os
from PySide6.QtCore import QThread, Signal


class Conversions(QThread):

    """
    A class that is responsible for carrying out the calculation of 
    permeability in the blood-brain barrier for a group of molecules in smi format
    """

    progress = Signal(int)
    finished = Signal()

    def __init__(self, ligand_path, folder_text):

        """
        Initializes the Conversions class with parameters for conversion.

        Parameters:
        -----------
        ligand_path : str
            Path to files in smi format for conversion
        folder_text : str
            The destination folder path where the converted files will be saved.
        """

        super().__init__()
        self._running = True
        self.maxim = 0
        self.contator = 0
        self.des_folder = folder_text + "/"
        self.ligand_path = ligand_path

    def run(self):

        # Create the destination folder
        os.makedirs(self.des_folder, exist_ok=True)

        # Read the smi file
        with open(self.ligand_path, "r") as f:
            lines = f.readlines()

        smile_list = []
        ID_list = []

        # Recover the smiles and ID's text and save it on their respective list
        for line in lines:
            if not self._running:
                break
            smile_list.append(line.split()[0])
            ID_list.append(line.split()[1])

        BBB_list = []
        Probab_list = []

        # Calculate the probability and the BBB and save it on their respective list
        for line in smile_list:
            if not self._running:
                break
            
            BBB_list.append(self.consensus_model(line)[0])
            Probab_list.append(self.consensus_model(line)[1])

            self.contator +=1
            self.progress.emit(self.contator)
        
        # Create a dictionary with the previous information
        My_dict = {"SMILES": lines,
                   "ID": ID_list,
                   "BBB": BBB_list,
                   "Probab": Probab_list}
        
        # Create a dataframe with the dictionary
        df = pd.DataFrame(My_dict)

        # Save the dictionary 
        df.to_csv(self.des_folder + "Results.csv", index = False)

        self.finished.emit()

    def consensus_model(self, smiles):
        try: 
            # Morgan Fingerprint
            mol = Chem.MolFromSmiles(smiles)
            fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            fingerprint_list = list(fingerprint)
            finger_array = np.array(fingerprint_list).reshape(1, -1)

            # DataFrame with name columns
            columns = [f'{i}' for i in range(finger_array.shape[1])]
            finger_df = pd.DataFrame(finger_array, columns=columns)

            # Import models
            SVM_model = load("codes/model_SVM.joblib")
            RF_model = load("codes/model_RF.joblib")
            XGB_model = load("codes/model_XGBoost.joblib")

            # Calculate predictions
            y1 = SVM_model.predict(finger_df)
            y1_prob = SVM_model.predict_proba(finger_df)

            y2 = RF_model.predict(finger_df)
            y2_prob = RF_model.predict_proba(finger_df)

            y3 = XGB_model.predict(finger_df)
            if y3 == 0:
                y3 = ["BBB+"]
            else:
                y3 = ["BBB-"]
            y3_prob = XGB_model.predict_proba(finger_df)

            consensus_probability = ((y1_prob[0][0]+y2_prob[0][0]+y3_prob[0][0])/3)
            clase = None

            # Return predictions
            if (y1 == ["BBB+"]) and (y2 == ["BBB+"]) and (y3 == ["BBB+"]):
                clase = "BBB+"
            else:
                clase = "BBB-"
            
            return clase, consensus_probability
        except Exception as e:
            return "Invalid Smiles", 0.0

    def Maximum(self, ruta: str):
        list_smiles = []

        with open(ruta, 'r') as archivo:
            lineas = archivo.readlines()

        for linea in lineas:
            list_smiles.append(linea.split()[0])

        self.maxim = len(list_smiles)
    
    def stop(self):
        self._running = False