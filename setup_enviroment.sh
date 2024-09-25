#!/bin/bash

# Check if conda is installed, if not, install it
if ! command -v conda &> /dev/null; then
    echo "Conda not found. Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3.sh
    bash Miniconda3.sh -b -p $HOME/miniconda3
    rm Miniconda3.sh
    export PATH="$HOME/miniconda3/bin:$PATH"
    echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
else
    echo "Conda is already installed."
fi

# Check if GUI_env environment exists, if not, create it
if conda info --envs | grep -q "GUI_env"; then
    echo "The conda environment 'GUI_env' already exists."
else
    echo "Creating the conda environment 'GUI_env'..."
    conda create --name GUI_env -y
fi

# Activate the environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate GUI_env

# Check if OpenBabel is installed, if not, install it
if ! conda list | grep -q "openbabel"; then
    echo "Installing OpenBabel..."
    conda install -c conda-forge openbabel=3.1.1 -y
else
    echo "OpenBabel is already installed."
fi

# Check if Smina is installed, if not, install it via Conda
if ! conda list | grep -q "smina"; then
    echo "Installing Smina via Conda..."
    conda install -c conda-forge smina=2020.12.10 -y
else
    echo "Smina is already installed."
fi

# Install the libraries required for the graphical interface.
echo "Installing the needed libraries..."
sudo apt update
sudo apt install -y libegl1 libopengl0 libgl1-mesa-glx

# Set permissions for project files
chmod -R +x "$(pwd)/project/lib/"
chmod +x "$(pwd)/project/Main_window"

# Reload ~/.bashrc to update paths
echo "Reloading ~/.bashrc..."
source ~/.bashrc

echo "The environment has been successfully configured. Done!"
