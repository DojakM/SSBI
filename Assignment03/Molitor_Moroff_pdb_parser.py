import glob

def read_residues(line):
    line = line[19:]
    return line.split(" ")

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
    return (type, length)

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
    return length

if __name__ == '__main__':
    helices = []
    sheets = []
    residues_list = []
    for path in glob.glob("Supplementary/**"):
        residues = []
        with open(path, 'r') as file:
            for line in file:
                if line.startswith("HELIX"):
                    helices.append(read_helix_entry(line))
                elif line.startswith("SHEET"):
                    sheets.append(read_sheet_entry(line))
                elif line.startswith("SEQRES"):
                    residues += read_residues(line)
        residues_list.append(residues)

    with open("output.txt", 'w') as file:
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
        str = "Secondary structure\tAbundance\nright-handed alpha helix\t{}\nright-handed 3_10 helices\t{}\nbeta-sheets" \
              "\t{}\n".format(alpha/total, hel_310/total, sheet_len/total)

        file.write(str)