import subprocess
import os
import glob
from Bio import PDB


class Conversions():
    def __init__(self):
        self.contator = 0
        self.maxim = 0
        self.chains = []
        self.selected_chains = []
        self.proteins = []
    
    def get_proteins(self, ini_folder):
        ini_folder = ini_folder + "/"
        pdbs_files = glob.glob(ini_folder + '*.pdb')

        for pdb in pdbs_files:
            self.proteins.append(os.path.basename(pdb))

    
    def Chains(self, ini_folder):
        ini_folder = ini_folder + "/"
        pdbs_files = glob.glob(ini_folder + '*.pdb')
        # Read and write pdb files
        parser = PDB.PDBParser(QUIET=True)

        for pdb in pdbs_files:

            # Read chains
            structure = parser.get_structure('protein', pdb)
            self.chains = []
            for model in structure:
                for chain in model:
                    self.chains.append(chain.id)
            
            yield

    def conversions(self, ini_folder, des_folder, save_ligands):
        ini_folder = ini_folder + "/"
        des_folder = des_folder + "/"
        os.makedirs(des_folder, exist_ok=True)
        pdbs_files = glob.glob(ini_folder + '*.pdb')

        # Read and write pdb files
        parser = PDB.PDBParser(QUIET=True)
        io = PDB.PDBIO()

        for pdb in pdbs_files:
            
            # Create the temporal file withouth the HETATM lines
            with open(pdb, 'r') as file:
                lines = file.readlines()
                num_line = [i for i, line in enumerate(lines) if "HETATM" in line]
                ligands = [line for line in lines if line.startswith("HETATM") and "HOH" not in line and "BHOH" not in line]
            
            # Filtered list
            filtered_list = [item for i, item in enumerate(lines) if i not in num_line]
            chains_ligands = {}
            
            # Read chains
            structure = parser.get_structure('protein', pdb)
            temporal_chains = []
            for model in structure:
                for chain in model:
                    temporal_chains.append(chain.id)

          
            # Save the ligands
            if save_ligands and ligands:
                ligands_folder = des_folder + "ligands_pdbqt/"
                os.makedirs(ligands_folder, exist_ok=True)

                # Dictionary for storing ligands by chain
                chains_ligands = {}

                for chain in temporal_chains:
                    # We filter the ligands corresponding to each chain without modifying the original list.
                    ligands_for_chain = [line for line in ligands if line.split()[4] == chain]
                    chains_ligands[chain] = ligands_for_chain

                # Save the ligands files in separated files and convert to pdbqt
                for key, value in chains_ligands.items():
                    if value:
                        pdb_file = ligands_folder + os.path.basename(pdb).split(".")[0] + "_ligand_" + key + ".pdb"
                        pdbqt_file = ligands_folder + os.path.basename(pdb).split(".")[0] + "_ligand_" + key + ".pdbqt"

                        # Save the ligands in .pdb format
                        with open(pdb_file, "w") as file:
                            for line in value:
                                file.write(line)

                        # Add End
                        # Save the ligands in .pdb format
                        with open(pdb_file, "a") as file:
                            file.write("END")

                        # Convert the .pdb file to .pdbqt using OpenBabel
                        result = subprocess.run(["./lib/obabel", "-i", "pdb", pdb_file, "-xr", "-opdbqt", '-O', pdbqt_file], check=True, capture_output=True, text=True)

                        # Remove the original .pdb file after conversion
                        os.remove(pdb_file)

            # Write the temporal files with only selected chains
            temp_pdb = pdb + "_temp"
            with open(temp_pdb, 'w') as temp_file:
                for line in filtered_list:
                    temp_file.write(line)

            structure_temp = parser.get_structure('protein', temp_pdb)
            io.set_structure(structure_temp)

            io.save(temp_pdb, select=ChainSelect(self.selected_chains))

            final_name = des_folder + os.path.basename(os.path.splitext(pdb)[0]) + ".pdbqt"
            subprocess.run(["./lib/obabel", "-i", "pdb", temp_pdb, "-xr", "-opdbqt", "-O", final_name], capture_output=True, text=True)
            
            os.remove(temp_pdb)

            self.contator += 1
            yield

    def Maximum(self, ini_folder):
        ini_folder = ini_folder + "/"
        pdbs_files = glob.glob(ini_folder + '*.pdb')
        self.maxim = len(pdbs_files)

class ChainSelect(PDB.Select):
    def __init__(self, chains_to_select):
        self.chains_to_select = chains_to_select

    def accept_chain(self, chain):
        # Aceptar solo las cadenas en la lista seleccionada
        return chain.id in self.chains_to_select
