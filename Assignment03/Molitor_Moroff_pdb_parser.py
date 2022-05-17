import glob
from functions import *

if __name__ == '__main__':
    for path in glob.glob("Supplementary/**"):
        with open(path, 'r') as file:
            for line in file:
                if line.startswith("HELIX"):
                    print(read_helix_entry(line))
                elif line.startswith("SHEET"):
                    print(read_sheet_entry(line))


