
import numpy


def contacts_intradomain(p, q, cont_matrix):
    """Calculates the number of interdomain connections.
    Parameters:
        p = number of residue at which the domain starts
        q = number of residue at which the domain ends
        cont_matrix = contact matrix of domain A and domain B, the location of close CA contacts are noted by 1"""
    # definition comes here for a domain starting at postion p and ending at position q.
    # The concats in the corresponding upper triangle matrix should be counted and returned.
    contacts_i = 0
    for x in range(p, q):
        for y in range(p, x):
            contacts_i += cont_matrix[x][y]
    return contacts_i


def contacts_interdomain(p1, p2, q1, q2, cont_matrix):
    """Calculates the number of interdomain connections between domain A and domain B.
    Parameters:
        p1 = starting residue of domain A
        p2 = starting residue of domain B
        q1 = ending residue of domain A
        q2 = ending residue of domain B
        cont_matrix = contact matrix of domain A and domain B, the location of close CA contacts are noted by 1"""
    # definition comes here. Contacts should be counted in a rectangular area where
    # domain A starts at p1 and ends at q1 and domain B from p2 to q2.
    contacts_e = 0
    for x in range(p1, q1):
        for y in range(p2, q2):
            contacts_e += cont_matrix[x][y]
    return contacts_e


def main():
    # storing coordinates in arrays ranging from 0 to
    # the number of residues
    coordx = []
    coordy = []
    coordz = []
    # residue numbers will also be stored as numbering in PDB files
    # rarely starts with 1, residues can be missing etc.
    pdb_resnum = []

    # cutoff value for CA-CA distances
    cutoff = 8
    # number of residues that have been read in
    nres = 0

    filepath = "c:/Marci/Egyetem/Medical Biotech MSc/3rd Semester/Integrated Structural Bioinformatics/3rd_Assignment_Domain_assignment/3izp.pdb"

    # Reading in the PDB file:
    pdbfile = open(filepath, "r")
    for line in pdbfile:
        # only ATOM lines will be considered
        if line.startswith("ATOM"):
            # obtaining the atom name and getting rid of spaces
            atom = line[12:15].replace(" ", "")
            # the chain identifier is 1 character
            chain = line[21:22]
            # processing only lines of CA atoms in chain A
            if atom == "CA" and chain == "E":
                # the append function adds a value to the array. The substring at given positions of the line
                # should be converted to a floating-point number.
                coordx.append(float(line[30:37]))
                coordy.append(float(line[38:45]))
                coordz.append(float(line[46:53]))
                # storing residue numbers, they are integers
                pdb_resnum.append(int(line[22:27]))
                # we have appended a value to each array, they should have equal size and the data
                # for the same residue should be accessible at the same index that corresponds to nres
                # (the first index is zero and we set nres to zero above). Now we increase nres to keep track.
                nres = nres + 1

    # Initialization: filling the contact matrix with zeroes
    contact_matrix = [[0 for h in range(nres)] for k in range(nres)]

    # Filling the matrix, calculating the distance of every residue pair i,j
    for i in range(1, nres):
        for j in range(0, i):
            distance = ((coordx[i] - coordx[j]) * (coordx[i] - coordx[j])) + ((coordy[i] - coordy[j]) * (coordy[i] - coordy[j])) + ((coordz[i] - coordz[j]) * (coordz[i] - coordz[j]))
            distance = numpy.sqrt(distance)

            if distance <= cutoff:
                contact_matrix[j][i] = 1
                contact_matrix[i][j] = contact_matrix[j][i]

    # Iterating over all possible cuts between domain A and domain B and calculating the number of intra
    # and interdomain connects.
    to_plot = []
    writefile = open("c:/Marci/Egyetem/Medical Biotech MSc/3rd Semester/Integrated Structural Bioinformatics/3rd_Assignment_Domain_assignment/3izp_output.txt", "w")
    for x in range(1, nres-1):
        intA = contacts_intradomain(0, x, contact_matrix)
        intB = contacts_intradomain(x+1, nres, contact_matrix)
        extAB = contacts_interdomain(0, x+1, x, nres, contact_matrix)
        d = (intA * intB) / (extAB**2)
        print("intA:", intA, " intB:", intB, "  extAB:", extAB, " Residue num: {} \t DOMAK: {}".format(pdb_resnum[x], d))
        to_plot.append(d)
        outp = "intA:" + str(intA) + " intB:" + str(intB) + " extAB:" + str(extAB) + " Residue num: {} \t DOMAK: {}".format(pdb_resnum[x], d) + "\n"
        writefile.write(outp)
        """print("A:", intA, " B:", intB, "  AB:", extAB, "{} \t {}".format(pdb_resnum[x], d))
        to_plot.append(d)
        outp = "A:" + str(intA) + " B:" + str(intB) + "  AB:" + str(extAB) + "{} \t {}".format(pdb_resnum[x], d) + "\n"
        writefile.write(outp)"""

main()
