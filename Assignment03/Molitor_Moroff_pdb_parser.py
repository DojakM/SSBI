import glob
import math
import matplotlib.pyplot as plt

def read_residues(line):
    line = line[20:]
    line = line.strip()
    values = line.split(" ")

    return values

def read_helix_entry(line):
    first_entry = line[22:25]
    last_entry = line[34:37]
    type = line[39:40]
    try:
        first_entry = int(first_entry)
        last_entry = int(last_entry)
        type = int(type)
    except:
        print("Incorrect read-in")
        quit()
    length = last_entry-first_entry
    return type, length, first_entry, last_entry

def read_sheet_entry(line):
    first_entry = line[23:26]
    last_entry = line[34:37]
    try:
        first_entry = int(first_entry)
        last_entry = int(last_entry)
    except:
        print("Incorrect read-in")
        quit()
    length = last_entry - first_entry
    return length, first_entry, last_entry

def read_atoms(line):
    if line[13] == "N" and line[14] == " ":
        coordinates_N = (float(line[31:38]),float(line[39:46]),float(line[47:53]))
        return coordinates_N, 1
    elif line[13] == "O" and line[14] ==" ":
        coordinates_O = (float(line[31:38]),float(line[39:46]),float(line[47:53]))
        return coordinates_O, 2
    else:
        return 0, 0

if __name__ == '__main__':
    helices = []
    sheets = []
    residues_list = []
    amino_acids = {"ALA":0, "ARG":1, "ASN":2, "ASP":3, "CYS":4, "GLU":5, "GLN":6,
                   "GLY":7, "HIS":8, "ILE":9, "LEU":10, "LYS":11, "MET":12, "PHE":13,
                   "PRO":14,  "SER":15,  "THR":16, "TRP":17,  "TYR":18,  "VAL":19}
    aa_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
                   "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                   "PRO",  "SER",  "THR", "TRP",  "TYR",  "VAL"]
    aa_propensity_1 = [0]*20
    aa_propensity_2 = [0] * 20
    length_distance = []
    #find once every file
    for path in glob.glob("Supplementary/**"):
        # many lists, probably a bit redundant
        read_in = False
        chain = "w"
        ns = []
        os = []
        offset = []
        hel_pla = []
        she_pla = []
        residues = []
        # for every file
        with open(path, 'r') as file:
            #for every line
            for line in file:
                # task a and b
                if line.startswith("HELIX"):
                    type, leng, first, last = read_helix_entry(line)
                    helices.append((type, leng)) #task a
                    if type == 1: #task b
                        hel_pla.append((first, last))
                elif line.startswith("SHEET"):
                    leng, first, last = read_sheet_entry(line)
                    sheets.append(leng) #task a
                    she_pla.append((first, last)) #task b
                elif line.startswith("SEQRES"):
                    residues += read_residues(line) #task b
                elif line.startswith("ATOM"):
                    if not read_in: #task b and c
                        offset.append(int(line[23:26]))
                        read_in = True
                        chain = line[21]
                    if chain != line[21]:
                        read_in = False
                        for i in range(4, len(ns)): #task c
                            n_coord = ns[i]
                            o_coord = os[i-4]
                            distance = math.sqrt((n_coord[0]-o_coord[0])**2+((n_coord[1]-o_coord[1])**2)+((n_coord[
                                                                                                               2]-o_coord[2])**2))
                            length_distance.append(distance)
                        ns = []
                        os = []
                    # task c
                    entry, mode = read_atoms(line)
                    if mode == 1:
                        ns.append(entry)
                    if mode == 2:
                        os.append(entry)
            for i in range(4, len(ns)): #task c
                n_coord = ns[i]
                o_coord = os[i - 4]
                distance = math.sqrt((n_coord[0] - o_coord[0]) ** 2 + ((n_coord[1] - o_coord[1]) ** 2) + ((n_coord[2] -o_coord[2]) ** 2))
                length_distance.append(distance)

        # Task 2b and c for helices run through and filter results
        residues_list.append(residues)
        counter = offset[0]
        run = 0
        heli_dist = []
        for i in hel_pla:
            # making sure to have the correct offset
            for j in range(i[0], i[1]):
                if len(offset)>run+1:
                    if counter + j >= offset[run+1]:
                        counter = offset[run+1]
                        run += 1
                if j > i[0]+4:
                    heli_dist.append(length_distance[j-counter-4])
                try:
                    aa_propensity_1[amino_acids.get(residues[j-counter])] = aa_propensity_1[amino_acids.get(residues[
                                                                                                               j-counter])]+ 1
                except:
                    continue
        # Task 2b only betasheet filter
        for i in she_pla:
            for j in range(i[0], i[1]):
                if len(offset) > run + 1:
                    if counter + j >= offset[run + 1]:
                        counter = offset[run + 1]
                        run += 1
                try:
                    aa_propensity_2[amino_acids.get(residues[j-counter])] = aa_propensity_2[amino_acids.get(residues[
                                                                                                               j-counter])]+1
                except:
                    continue

    # Task 2b Output txt file with the propensity of one amino acid being in any givin structure
    with open("output2.txt", 'w') as file:
        total1 = sum(aa_propensity_1)
        total2 = sum(aa_propensity_2)
        string = "Amino acids\talpha-helix\tbeta-sheet\n"

        for i in range(0, 20):
            percent1 = aa_propensity_1[i]/total1
            percent2 = aa_propensity_2[i]/total2
            string += "{}\t{}\t{}\n".format(aa_list[i], str(percent1), str(percent2))
        file.write(string)

    #Task 2a write a txt file with the length of all combined structures of one type/total of structure lengths
    with open("output.txt ", 'w') as file:
        sheet_len = sum(sheets)
        helices_len = 0
        alpha = 0
        hel_310 = 0
        for i in helices:
            helices_len += i[1]
            if i[0] == 1:
                alpha += i[1]
            if i[0] == 5:
                hel_310 += i[1]
        total = 0
        for residues in residues_list:
            total += len(residues)
        string = "Secondary structure\tAbundance\nright-handed alpha helix\t{}\nright-handed 3_10 helices\t{}\nbeta-sheets\t{}\n".format(alpha/total, hel_310/total, sheet_len/total)
        file.write(string)

    # Task 2c plotting of the histograms
    fig, ax = plt.subplots(1,2)
    ax[0].hist(heli_dist, bins=20, density=True, color="red")
    ax[1].hist(length_distance, bins=200, density=True, color="green")
    plt.title("Comparison of molecule distances in (left) helices and (right) all structures")
    plt.savefig("Assignment03.png")
