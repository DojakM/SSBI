import glob
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
    if line[18] == "N":
        return None
    if line[18] == "O":
        return None


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
    for path in glob.glob("Supplementary/**"):
        read_in = False
        chain = "w"
        offset = []
        hel_pla = []
        she_pla = []
        residues = []
        with open(path, 'r') as file:
            for line in file:
                if line.startswith("HELIX"):
                    type, leng, first, last = read_helix_entry(line)
                    helices.append((type, leng))
                    if type == 1:
                        hel_pla.append((first, last))
                elif line.startswith("SHEET"):
                    leng, first, last = read_sheet_entry(line)
                    sheets.append(leng)
                    she_pla.append((first, last))
                elif line.startswith("SEQRES"):
                    residues += read_residues(line)
                elif line.startswith("ATOM"):
                    if not read_in:
                        offset.append(int(line[23:26]))
                        read_in = True
                        chain = line[22]
                    if chain != line[22]:
                        read_in = False
                    read_atoms(line)
        residues_list.append(residues)
        counter = offset[0]
        run = 0
        print(len(offset))
        for i in hel_pla:
            for j in range(i[0], i[1]):
                if len(offset)>run+1:
                    if counter + j >= offset[run+1]:
                        counter = offset[run+1]
                        run += 1
                try:
                    aa_propensity_1[amino_acids.get(residues[j-counter])] = aa_propensity_1[amino_acids.get(residues[
                                                                                                               j-counter])]+ 1
                except:
                    continue

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

        with open("output2.txt", 'w') as file:
            total1 = sum(aa_propensity_1)
            total2 = sum(aa_propensity_2)
            string = "Amino acids\talpha-helix\tbeta-sheet\n"

            for i in range(0, 20):
                percent1 = aa_propensity_1[i]/total1
                percent2 = aa_propensity_2[i]/total2
                string += "{}\t{}\t{}\n".format(aa_list[i], str(percent1), str(percent2))
            file.write(string)

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