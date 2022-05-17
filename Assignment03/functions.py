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