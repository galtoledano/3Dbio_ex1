# 3D-Bio, Ex1
# by yahel yedidya and gal toledano

import Bio.PDB as pdb
import sys
from Bio import pairwise2


# BONUS 9.2 IMPLEMENTED
def align_structs(id1, chain1, id2, chain2):
    """
    the main function. gets the ids and the chain's names and finds the alignment with the best RMSD.
    prints the best RMSD, and saving the alignments file in cif format
    :param id1: the first file id
    :param chain1: the first protein's chain
    :param id2: the second file id
    :param chain2: the second protein's chain
    """
    # generating the relevant data
    lst = pdb.PDBList()
    protein1 = lst.retrieve_pdb_file(id1)
    protein2 = lst.retrieve_pdb_file(id2)
    parser = pdb.MMCIFParser()
    struct1 = parser.get_structure("p1", protein1)
    struct2 = parser.get_structure("p2", protein2)

    # creating a lists of CA atoms to align
    atoms1 = create_atoms_list(struct1, chain1)
    atoms2 = create_atoms_list(struct2, chain2)
    if len(atoms1) != len(atoms2):
        atoms1, atoms2 = bonus_9_2(chain1, chain2, struct1, struct2)

    # making the align
    super_imposer = pdb.Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(struct2[0].get_atoms())
    print(super_imposer.rms)

    # saving the aligned structure to files
    saving_file(id1, struct1)
    saving_file(id2, struct2)


def bonus_9_2(chain1, chain2, struct1, struct2):
    """
    handling the case that the number of atoms is different by align the amino-acids and remove the unnecessary atoms.
    :param chain1: the first protein's chain
    :param chain2: the second protein's chain
    :param struct1: the first protein
    :param struct2: the second protein
    :return: filtered lists of atoms to align
    """
    ppbils = pdb.CaPPBuilder()
    peptide1, peptide2 = [], []
    filter_peptide_by_chain(chain1, peptide1, ppbils, struct1)
    filter_peptide_by_chain(chain2, peptide2, ppbils, struct2)

    # converting list to peptide
    peptide1 = pdb.Polypeptide.Polypeptide(peptide1)
    peptide2 = pdb.Polypeptide.Polypeptide(peptide2)

    # converting peptide to sequence and CA list
    seq1 = peptide1.get_sequence()
    atoms1 = peptide1.get_ca_list()
    seq2 = peptide2.get_sequence()
    atoms2 = peptide2.get_ca_list()

    # align the sequences
    alignments = pairwise2.align.globalxx(seq1, seq2)

    # filter the atoms lists
    ignore_inx1 = [i for i in range(len(alignments[0][0])) if alignments[0][0][i] == "-"]
    ignore_inx2 = [i for i in range(len(alignments[0][1])) if alignments[0][1][i] == "-"]
    atoms2 = [atoms2[i] for i in range(len(atoms2)) if i not in ignore_inx1]
    atoms1 = [atoms1[i] for i in range(len(atoms1)) if i not in ignore_inx2]
    return atoms1, atoms2


def filter_peptide_by_chain(chain, peptide, ppbils, struct):
    """
    this function filters the peptide atoms by given chain
    :param chain: the given chain
    :param peptide: the filtered list
    :param ppbils: the peptide to filter
    :param struct: the structure
    """
    for pp in ppbils.build_peptides(struct):
        for r in pp:
            if r.get_parent().get_id() == chain:
                peptide.append(r)


def saving_file(id, struct):
    """
    saving the file in cif format
    :param id: the file's id
    :param struct: the structure we want to save
    """
    io = pdb.MMCIFIO()
    io.set_structure(struct)
    io.save("{0}.cif".format(id))


def create_atoms_list(struct, chain_name):
    """
    creates a list with CA atoms to align from the given structure at the given chain
    :param struct: the structure
    :param chain_name: the chain we want to align
    :return: the atoms as list
    """
    lst = []
    for model in struct:
        for chain in model:
            if chain.get_id() == chain_name:
                lst += [atom for atom in chain.get_atoms() if atom.get_id() == 'CA']
    return lst


if __name__ == '__main__':
    # prints Error in case of invalid number of arguments
    if len(sys.argv) != 5:
        print("Error: invalid number of arguments")
    align_structs(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])



