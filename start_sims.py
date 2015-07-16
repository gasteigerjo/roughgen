import pickle
import glob
import shutil
import os
import subprocess

if __name__ == '__main__':
    # Get directory names
    dirs_file = open("dirs.txt",'r')
    dirs = pickle.load(dirs_file)
    dirs_file.close()

    # For each mesh created
    oldDir = os.getcwd()
    for directory in dirs:
        # Copy the necessary files from setup
        if os.path.exists(directory):
            shutil.rmtree(directory)
        shutil.copytree("{0}/setup/".format(oldDir), directory, symlinks=True)

        # Copy the mesh
        shutil.copytree("{0}_temp".format(directory), "{0}/mesh".format(directory), symlinks=True)
        shutil.rmtree("{0}_temp".format(directory))

        # Run the simulation
        os.chdir(directory)
        subprocess.call("mpiexec.hydra -n $mpi_ranks ./SeisSol_release_generatedKernels_dsnb_hybrid_none_9_5 PARAMETERS.par", shell=True)
    os.chdir(oldDir)
