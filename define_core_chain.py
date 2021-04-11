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
