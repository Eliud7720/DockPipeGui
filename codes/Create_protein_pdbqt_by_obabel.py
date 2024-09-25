import subprocess
import os
import glob
from Bio import PDB
from PySide6.QtCore import QThread, Signal


class Conversions(QThread):
    progress = Signal(int)
    finished = Signal()

    """
    A function responsible for generating files in pdbqt format from a 
    folder with proteins files on pdb format using openbabel
    """

    def __init__(self, ini_folder, des_folder, save_ligands):

        """
        Initializes the Conversions class with parameters for conversion.

        Parameters:
        -----------
        ini_folder : str
            Path to the folder with proteins in pdb format for the conversion to pdbqt
        des_folder : str
            The destination folder path where the converted files will be saved.
        save_ligands : int
            A number that indicates if save the ligands 
        """

        super().__init__()
        self.contator = 0
        self.maxim = 0
        self.chains = []
        self.selected_chains = []
        self.proteins = []
        self._running = True
        self.ini_folder = ini_folder + "/"
        self.des_folder = des_folder + "/"
        self.save_ligands = save_ligands
    
    def get_proteins(self, ini_folder):

        """
        A method that determines the proteins present in the folder
        """

        ini_folder = ini_folder + "/"
        pdbs_files = glob.glob(ini_folder + '*.pdb')

        for pdb in pdbs_files:
            self.proteins.append(os.path.basename(pdb))

    def Chains(self, ini_folder):

        """
        A method that determines the chains of each protein
        """

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
            
            yield self.chains

    def run(self):

        # Make the destination folder
        os.makedirs(self.des_folder, exist_ok=True)
        pdbs_files = glob.glob(self.ini_folder + '*.pdb')

        # Read and write pdb files
        parser = PDB.PDBParser(QUIET=True)
        io = PDB.PDBIO()

        for pdb in pdbs_files:

            if not self._running:
                break
            
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
            if self.save_ligands and ligands:
                ligands_folder = self.des_folder + "ligands_pdbqt/"
                os.makedirs(ligands_folder, exist_ok=True)

                # Dictionary for storing ligands by chain
                chains_ligands = {}

                for line in ligands:
                    if line.split()[4] not in chains_ligands.keys():
                        chains_ligands[line.split()[4]] = []

                    chains_ligands[line.split()[4]].append(line)

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
                        result = subprocess.run(["obabel", "-i", "pdb", pdb_file, "-xr", "-opdbqt", '-O', pdbqt_file, "-p", "7.0", "--partialcharge", "gasteiger"], check=True, capture_output=True, text=True)

                        with open(pdbqt_file, 'r') as file:
                            lines = file.readlines()
                            lines.insert(0, "ROOT\n")
                            lines.insert(-1, "ENDROOT\n")
                            lines.insert(-1, "TORSDOF 0\n")
                        
                        with open(pdbqt_file, 'w') as file:
                            file.writelines(lines)

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

            final_name = self.des_folder + os.path.basename(os.path.splitext(pdb)[0]) + ".pdbqt"

            # Run the command to convert the protein with the selected chains to pdbqt
            result = subprocess.run(["obabel", "-i", "pdb", temp_pdb, "-xr", "-opdbqt", "-O", final_name, "--partialcharge", "gasteiger"], capture_output=True, text=True)
            
            os.remove(temp_pdb)

            self.contator += 1
            self.progress.emit(self.contator)
        
        self.finished.emit()

    def Maximum(self, ini_folder):
        ini_folder = ini_folder + "/"
        pdbs_files = glob.glob(ini_folder + '*.pdb')
        self.maxim = len(pdbs_files)
    
    def stop(self):
        self._running = False

class ChainSelect(PDB.Select):
    def __init__(self, chains_to_select):
        self.chains_to_select = chains_to_select

    def accept_chain(self, chain):
        return chain.id in self.chains_to_select
