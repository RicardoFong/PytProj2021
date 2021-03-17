from Bio.PDB import *
from Bio import pairwise2 # This library is for applying pairwise comparison of sequences
# from Bio.SubsMat import MatrixInfo as matlist # Used for applying subtitution matrix in sequence comparison
from Bio.pairwise2 import format_alignment # DEBUGING Used for printing the alignment. DEBUGING purposes
from math import log # For use in defining the gap function for sequence comparison parameters
import sys
import os
import argparse
import re

###################################################################################################################

parser = argparse.ArgumentParser(description = """This program analyzes pdb binary chain interactions provided,
                                                  calculate and reconstructs the polipeptide complexrepresented
                                                  by the provided files.""")

parser.add_argument('-i', '--input',
                                        dest = "inPath",
                                        action = "store",
                                        default = os.getcwd(), # Default is current directory
                                        help = """Provide the complete path to the pdb files you
                                                  want to analyze after -i flag. If no path is provided,
                                                  the program will assume that current directory is
                                                  going to be searched for pdb files to analize.""")


parser.add_argument('-o', '--output',
                                        dest = "outfile",
                                        action = "store",
                                        default = None,
                                        help = """The output file will be named after the name provided after
                                                  this flag. This option is mandatory, if no output filename
                                                  is provided, no analysis will be run and this message will
                                                  be displayed. """)

parser.add_argument('-v', '--verbose',
                                        dest = "verbose",
                                        action = "store_true",
                                        default = False,
                                        help = """This option prints the execution of the program as it is being
                                                  performed. It is adviseable to run the analysis with this option
                                                  in case any exception is raised.""")

args = parser.parse_args()

###################################################################################################################

def get_files(input):
    """ This method reads the directory that is provided in the command line argument.
        Default directory to read is the current directory. It searches for files with
        .pdb extension and returns them."""
    path = input
    extension = re.compile('.pdb$') # create a regular expression object refering to .pdb file extension

    # list comprehension, loops through all files in "input" (input directory)
    pdb_files = [file for file in os.listdir(path) if extension.search(file) is not None]
    # changes working directory to input directory
    os.chdir(path)

    # returns the names of files with .pdb extension
    return pdb_files

def get_file_prefix(pdb_files):
        """Tests for match of regex containing .pdb extension, trims it and returns
           the file prefix."""

        p = re.compile('(.*).pdb$') # compile regular expression matching any file with .pdb extension
        m = p.search(pdb_files) # searches for match in the file names provided

        return m.group(1) # returns just the prefix of files matching regex

###################################################################################################################
# def get_seq(chain):
    # """ gets sequence from chain provided. Uses CaPPBuilder from Bio.pdb, the distance between the
    #     alpha carbons can be changed inserting radius = integer or float in CaPPBuilder."""


#     ppb = CaPPBuilder() # set parser for alpha carbon
#     parser = PDBParser(PERMISSIVE = 1, QUIET=True) # set file parser
#
#     for pp in ppb.build_peptides(chain): # looping through alpha carbons of chain
#         seq = pp.get_sequence() # getting sequence, aplha carbons
#         yield seq
###################################################################################################################

def compare_seqs(seq1, seq2):
    """This method aligns two given sequences, calculates percentage identity and
       tests if this identity is equal or higher than 95% or lower. If the given
       sequences show 95% identity or more, method returns True. otherwise it returns
       False."""

    # definition of BLAST percentage identity is taken from:
    # Heng Li, 2021, at: https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity

    # value 1 is for counting identical matches. Therefore the alignment score is
    # equal to counting identical matches. percentage identity is then calculated
    # as the ratio between the score and the length of the alignment. this percentage
    # identity is the BLAST definition.
    align = pairwise2.align.globalmx(seq1, seq2, 1, 0) # setting the alignment in align object
    algn_len = len(align[0][1]) # getting the length of the alignment with gaps
    match_score = align[0][2] # getting the score of alignment

    identity = match_score/algn_len # percentage identity is ratio of identical matches and alignment length

    # evaluating if pair of sequences are homodimers or heterodimers
    # True = homodimer. False = heterodimer
    if identity >= .95:
        return True
    else:
        return False

def get_structures(pdb_files):
    """ This method parses pdb files. It loops over structures and chains of the
        pdb files and returns a series of lists, each one containing the polipeptide
        sequence of the chains contained in the pdb files. Input is an iterator
        comprised of pdb files. if any of the chains fail to contain alpha carbon
        interactions, the method will represent it by a single 'X'."""

    # setting parser of alpha carbons and pdb file parser in objects
    ppb = CaPPBuilder()
    parser = PDBParser(PERMISSIVE = 1, QUIET=True)

    # looping through files provided in input
    for file in pdb_files:

        # sequence variables to be used for getting polipeptide sequences from chains
        seq = []
        seq1 = []
        sequences = []
        sequences1 = []
        chains = []

        id = get_file_prefix(file) # id is the prefix of input files
        structure = parser.get_structure(id, file) # set pdb file structure in structure object
        models = structure[0] # setting models in object

        # looping through models
        for chain in models:

            # filling list of chains
            chains.append(chain)

        # looping through first chain of all pdb files
        for pp in ppb.build_peptides(chains[0]):

            extract = pp.get_sequence() # extracting sequence of alpha carbons
            seq.append(str(extract)) # filling seq list with sequences as strings
            sequences = ["".join(seq)] # if chains yield more than one seq fragment join them into one string

        # looping through second chain of all pdb files
        for pp in ppb.build_peptides(chains[1]):

            extract = pp.get_sequence()
            seq1.append(str(extract))
            sequences1 = ["".join(seq1)]

        # testing if both sequences are empty (in which case they don't belong to a polipeptide)
        if len(sequences1) == 0:
            sequences1 = ["X"]
        if len(sequences) == 0:
            sequences = ["X"]

        merge = zip(sequences,sequences1) # setting an object with both sets of sequences, zipped
        seq_tuples = list(merge) # creating an object containing all  tuples of pairs of sequences

        # yields a generator object with all the lists of tuples, containing pairs of sequences
        yield seq_tuples

###################################################################################################################

if __name__=="__main__":

    # creating object with all the files in the directory provided
    files = get_files(args.inPath)
    # getting structure of files, while getting their pdb structure
    for file in get_structures(files):
        # unraveling tuples of sequences from file structures
        for iter in file:

            if iter[1] == 'X' and iter[0] == 'X':
                pass # if both sequences are empty, do nothing
            else:
                u = compare_seqs(iter[0], iter[1]) # comparing sequences in tuple
                print(u)
###################################################################################################################
########## Not used for project. But useful for calculating protein evolutionary distance #########################
###################################################################################################################

# def gap_opening_function(x, y): # x is gap position in seq, y is gap length
# """ This function handles the penalties of gap openings and extensions in
# compare_seqsBLOSUM method. change the values on the return funtions for changing
# the penalty values. """
#
# if y == 0: # No gap
# return 0
# elif y == 1: #Gap open penalty
# return -2
# return - (2 + y/4.0 + log(y)/2.0)
#
# def gap_extension_function(x, y): # x is gap position in seq, y is gap length
# """ This function handles the penalties of gap openings and extensions in
# compare_seqsBLOSUM method. change the values on the return funtions for changing
# the penalty values. """
#
# if y == 0: # No gap
# return 0
# elif y == 1: #Gap open penalty
# return -2
# return - (2 + y/4.0 + log(y)/2.0)
#
# def compare_seqsBLOSUM50(seq1, seq2):
#     """This method returns the score of a pairwise sequence comparison. Expects
#     two sequences. The matrix substitution is blosum50, and penalizes gap openings
#     and extensions with values set in gap_function. default penalties are :
#
#     gap opening: -2 # I am exploring which penalties are best for gap opening and extension.
#     gap extension: -2 """
#
#     matrix = matlist.blosum50  # Bio.Align.substitution_matrices   -   I need to set the substitution matrix with this module instead
#     alignment = pairwise2.align.globaldc(seq1, seq2, matrix, gap_opening_function, gap_extension_function)
#     score = alignment[0][2]
###################################################################################################################
###################################################################################################################
