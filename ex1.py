import Bio.PDB as pdb
import sys
from Bio import pairwise2
import Bio.Seq as seq


def q9_1(id1, chain1, id2, chain2):
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

    # ppbils = pdb.CaPPBuilder()
    # pp_lst1, pp_lst2 = seq.Seq(""), seq.Seq("")
    # for pp in ppbils.build_peptides(struct1):
    #     pp_lst1 += pp.get_sequence()
    # for pp in ppbils.build_peptides(struct2):
    #     pp_lst2 += pp.get_sequence()
    # alignments = pairwise2.align.globalxx(pp_lst1, pp_lst2)
    # print(alignments)
    # ignore_inx1 = [i for i in range(len(alignments[0][0])) if alignments[0][0][i] == "-"]
    # ignore_inx2 = [i for i in range(len(alignments[0][1])) if alignments[0][1][i] == "-"]

    # creating atoms to align lists
    atoms1 = create_atoms_list(struct1, chain1, [])
    atoms2 = create_atoms_list(struct2, chain2, [])
    # if len(atoms1) != len(atoms2):
    #     atoms2 = [atoms2[i] for i in range(len(atoms2)) if i not in ignore_inx1]
    #     atoms1 = [atoms1[i] for i in range(len(atoms1)) if i not in ignore_inx2]

    # making the align
    super_imposer = pdb.Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(struct2[0].get_atoms())
    print(super_imposer.rms)

    # saving to file
    saving_file(id1, struct1)
    saving_file(id2, struct2)


def saving_file(id, struct):
    """
    saving the file in cif format
    :param id: the file's id
    :param struct: the structure we want to save
    """
    io = pdb.MMCIFIO()
    io.set_structure(struct)
    io.save("{0}.cif".format(id))


def create_atoms_list(struct, chain_name, ignore_inx):
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
                lst += [atom for atom in chain.get_atoms() if atom.get_name() == 'CA']
                # for i in range(len(atoms)):
                #     if atoms[i].get_name() == 'CA' and i not in ignore_inx:
                #         lst.append(atoms[i])
    return lst


if __name__ == '__main__':
    # prints Error in case of invalid number of arguments
    if len(sys.argv) != 5:
        print("Error: invalid number of arguments")
    q9_1(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])



