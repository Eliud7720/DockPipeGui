import subprocess
import os
import glob
import MDAnalysis as mda
import prolif as plf
from rdkit import Chem
import numpy as np
import matplotlib.pyplot as plt

class Conversions():
    def __init__(self):
        self.maxim = 0
        self.contator = 0

    def conversions(self, protein_file, pdbqt_files, folder_text):
        des_folder = folder_text + "/"
        ligands = pdbqt_files + "/"
        os.makedirs(des_folder, exist_ok=True)
        ligands_files = glob.glob(ligands + '*.pdbqt')
        ligands_files = sorted(ligands_files)

        # add H to the protein
        result = subprocess.run(["./lib/obabel", protein_file, "-O", des_folder + os.path.basename(protein_file), "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Prepare the protein
        u = mda.Universe(des_folder + os.path.basename(protein_file), topology_format="PDB", format="PDB")
        protein_mol = plf.Molecule.from_mda(u)

        # Convert PDBQT files to smiles
        smiles_list = []
        for ligand in ligands_files:
            result = subprocess.run(["./lib/obabel", ligand, "-O", des_folder + "ligand.pdb"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            mol = Chem.MolFromPDBFile(des_folder + "ligand.pdb")
            smiles = Chem.MolToSmiles(mol)
            smiles_list.append(smiles)

            os.remove(des_folder + "ligand.pdb")
        
        template_list = []    
        for i, ligand in enumerate(ligands_files):
            template = Chem.MolFromSmiles(smiles_list[i])
            template_list.append(template)
        
        pose_iterable = plf.pdbqt_supplier(ligands_files, template_list)
        fp = plf.Fingerprint()

        fp.run_from_iterable(pose_iterable, protein_mol)
        df = fp.to_dataframe(index_col="Molecule")
        df.to_csv(des_folder + "Interactions results.csv")

        # Save fig
        fp.plot_barcode(xlabel="Ligand")
        plt.savefig(des_folder + "Interactions_graphic.png", format='png', dpi=300)
        plt.close()

        for i in range(0, len(pose_iterable)):
            fig = fp.plot_lignetwork(pose_iterable[i])
            with open(des_folder + f"Ligand_{i}.html", "w") as f:
                f.write(fp.plot_lignetwork(pose_iterable[0]).data)
    
        


    def Maximum(self, ini_folder):
        ini_folder = ini_folder + "/"
        pdbs_files = glob.glob(ini_folder + '*.pdbqt')
        self.maxim = len(pdbs_files)