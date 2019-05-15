import utilities as util

def flt(vcf_file, output):
    vcf, comments = util.readVcf(vcf_file)
    for line in vcf:
        if line[0] == "#":
            continue
        a = line.split()
        tmp = a[9].split(':')
        mother = tmp[0]
        tmp = a[10].split(':')
        father = tmp[0]
        tmp = a[11].split(':')
        son = tmp[0]
        is_filtered = False
        if mother == "0/0" and father == "0/0" and son == "1/1":
            is_filtered = True
        elif mother == "0/0" and father == "0/1" and son == "1/1":
            is_filtered = True
        elif mother == "0/1" and father == "0/0" and son == "1/1":
            is_filtered = True
        elif mother == "1/1" and father == "1/1" and son == "0/0":
            is_filtered = True
        elif mother == "0/1" and father == "1/1" and son == "0/0":
            is_filtered = True
        elif mother == "1/1" and father == "0/1" and son == "0/0":
            is_filtered = True
        elif mother == "0/0" and father == "1/1" and son == "0/0":
            is_filtered = True
        elif mother == "0/0" and father == "1/1" and son == "1/1":
            is_filtered = True
        elif mother == "1/1" and father == "0/0" and son == "1/1":
            is_filtered = True
        elif mother == "1/1" and father == "0/0" and son == "0/0":
            is_filtered = True

        if is_filtered == False:
            output.write("%s"%line)
    print("Done...")
    output.close()
