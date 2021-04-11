from statistics import mode # used for getting most comon element in list, for define_core_chain method
from Bio.PDB import *
from Bio import pairwise2 # This library is for applying pairwise comparison of sequences
# from Bio.SubsMat import MatrixInfo as matlist # Used for applying subtitution matrix in sequence comparison
from Bio.pairwise2 import format_alignment # DEBUGING Used for printing the alignment. DEBUGING purposes
from math import log # For use in defining the gap function for sequence comparison parameters
import sys
import os
import argparse
import re

######################################################################################################################################


####### I saw this first function is already in main program #################


####### define_core_chain is the bit of code needed to complete program ######


def chain_getter(pdb_files):
    """ expects file list as object from function get_files, outputs a generator containing a list of all chain occurrences """

    parser = PDBParser(PERMISSIVE = 1, QUIET=True)

    chains = []

    # looping through files provided in input
    for file in pdb_files:

        id = get_file_prefix(file)
        structure = parser.get_structure(id, file) # set pdb file structure in structure object

        for model in structure:
            for chain in model:
                yield chain

def define_core_chain(chain_iterator):
    """ expects a generator object containing the complete list of chains in a folder. returns the chain object in bio.pdb format with most occurrences """
    chains = []
    for element in chain_iterator:

        chains.append(element)
    return mode(chains)



############# grouping files of a same complex in a value list with complex name as key ##############

def complex_division(pdb_files):
    """ Parses a list object with names of files, returns a dictionary with the complex name as keys and correspondent files as values. or protein name as keys and files as values, depending on file name structure."""

    cmplx_name = []
    prot_name = []
    regex = re.compile("(\w+)\.(\w+)\.(\w+)_(\w)_(\w+)\.pdb")
    names_unique = []

    complexes = dict()

    for file in pdb_files:

        if regex.match(file):

            file_name_cmplx = get_file_prefix(file)
            cmplx = file_name_cmplx.split(".")[0]
            cmplx_name.append(cmplx)
            names_cmplx = set(cmplx_name)
            names_unique_cmplx = list(names_cmplx)
            cmplx_size = len(names_unique_cmplx)
        else:

            file_name_prot = get_file_prefix(file)
            protein = file_name_prot.split(".")[0]
            prot_name.append(protein)
            names_prot = set(prot_name)
            names_unique_prot = list(names_prot)
            prot_size = len(names_unique_prot)

    if regex.match(file):

        i = 0
        while i < cmplx_size:

            complexes[names_unique_cmplx[i]] = []

            i+=1

        for file in pdb_files:
            prefix = file.split(".")[0]
            if prefix in file:
                complexes.setdefault(prefix,[]).append(file)

    else:

        for file in pdb_files:
            prefix = file.split("_")[0]
            if prefix in file:
                complexes.setdefault(prefix,[]).append(file)


    return complexes




############ Sintaxis para realizar función principal por cada complex ##########################


if __name__=="__main__":

    # creating object with all the files in the directory provided


    files = get_files(args.inPath) # Ésta es la variable de archivos en mi script, solo hay que cambiarl por la pertinente.

    # calls function that builds a dictionary with name of complex as keys and file names as values
    files_dict = complex_division(files)

    # getting the identifier of the dictionary keys
    keys = list(files_dict.keys())
    complex = len(keys)


    i = 0 # initializing counter

    while i < complex: # make sure that program loops over each complex in the folder


        # all the actions for files executed in this section
        # will be executed for each complex in the folder separately, one by one

        #### example, see if it works: UNCOMMENT TO TEST

        # chains = chain_getter(files_dict[keys[i]])

        # core = define_core_chain(chains)

        # print(core)
