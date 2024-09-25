# Filter warning of prolif
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import subprocess
import os
import glob
import MDAnalysis as mda
import prolif as plf
from rdkit import Chem
import matplotlib.pyplot as plt
from PySide6.QtCore import QThread, Signal

class Conversions(QThread):

    """
    A class responsible for generating interactions files from a 
    folder with ligands files on pdbqt format using prolif
    """

    progress = Signal(int)
    finished = Signal()

    def __init__(self, protein_file, pdbqt_files, folder_text):

        """
        Initializes the Conversions class with parameters for conversion.

        Parameters:
        -----------
        protein_file : str
            Path to the file of the protein in pdb format
        pdbqt_files : str
            Path to the ligands folder in pdbqt format
        folder_text: str
            The destination folder path where the converted files will be saved.
        """

        super().__init__()
        self._running = True
        self.maxim = 0
        self.contator = 0
        self.protein_file = protein_file
        self.des_folder = folder_text + "/"
        self.ligands = pdbqt_files + "/"
        

    def run(self):

        # Make the destination folder
        os.makedirs(self.des_folder, exist_ok=True)

        # Search the pdbqt files 
        ligands_files = glob.glob(self.ligands + '*.pdbqt')
        ligands_files = sorted(ligands_files)

        # add H to the protein
        result = subprocess.run(["obabel", self.protein_file, "-O", self.des_folder + os.path.basename(self.protein_file), "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Prepare the protein
        u = mda.Universe(self.des_folder + os.path.basename(self.protein_file), topology_format="PDB", format="PDB")
        protein_mol = plf.Molecule.from_mda(u)

        # Convert PDBQT files to smiles using openbabel for each molecule
        smiles_list = []
        for ligand in ligands_files:

            if not self._running:
                break
            
            # Read the pdbqt file
            with open(ligand, 'r') as f2:
                lines = f2.readlines()
            
            # If the pdbqt file has the smiles, take that one, otherwise, calculate it
            if len(lines) > 3 and len(lines[3].split()) > 2:
                smiles = lines[3].split()[2]
            else:
                subprocess.run(["obabel", ligand, "-O", self.des_folder + "temp.pdb"], stdout=subprocess.DEVNULL,  stderr=subprocess.DEVNULL)
                pdb = Chem.MolFromPDBFile(self.des_folder + "temp.pdb")
                smiles = Chem.MolToSmiles(pdb)

            # Append each smiles to the smile list
            smiles_list.append(smiles)
            os.remove(self.des_folder + "temp.pdb")
            self.contator+=1
        
        # Create the template list for the run_from iterable (needs the mol object)
        template_list = []    
        for i, ligand in enumerate(ligands_files):
            if not self._running:
                break

            template = Chem.MolFromSmiles(smiles_list[i])
            template_list.append(template)
        
        # Create the pose iterable for each molecule
        pose_iterable = plf.pdbqt_supplier(ligands_files, template_list)
        fp = plf.Fingerprint()

        # Run the fp of prolif
        fp.run_from_iterable(pose_iterable, protein_mol)
        df = fp.to_dataframe(index_col="Molecule")
        df.to_csv(self.des_folder + "Interactions results.csv")

        # Save the fig
        fp.plot_barcode(xlabel="Ligand")
        plt.savefig(self.des_folder + "Interactions_graphic.png", format='png', dpi=300)
        plt.close()

        # Create the prolif html file for each file
        for i in range(0, len(pose_iterable)):
            if not self._running:
                break
            
            fig = fp.plot_lignetwork(pose_iterable[i])
            print(fig)
            print(pose_iterable[i])

            with open(self.des_folder + f"Ligand_{i}.html", "w") as f:
                f.write(fp.plot_lignetwork(pose_iterable[i]).data)

            self.contator +=1 
            self.progress.emit(self.contator)

        os.remove(self.des_folder + os.path.basename(self.protein_file))

        self.finished.emit()
        
    def Maximum(self, ini_folder):
        ini_folder = ini_folder + "/"
        pdbs_files = glob.glob(ini_folder + '*.pdbqt')
        self.maxim = len(pdbs_files) * 2

    def stop(self):
        self._running = False