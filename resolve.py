import utilities as util

def resolve(genome, vcf_file, vcf_new, mei_annot_file, is_precise):

    if is_precise:
        print("Sequence resolving the VCF file for only the precise calls")
    else:
        print("Sequence resolving the VCF file for all the calls")

    vcf, comments = util.readVcf(vcf_file)
    chr_name_old = "1"
    chr_change=False
    fa = util.readGenome(genome,chr_name_old)

    for line in comments:
        vcf_new.write("%s"%line)

    meis = dict()

    for line in vcf:
        a = line.split()
        if "INS:ME" in a[4]:
            tmp = a[4].split(':')
            mei_annot=tmp[2][:-1]
            meis.update({mei_annot:""})

    print("Processing chromosome",chr_name_old)
    for line in vcf:
        if ("IMPRECISE" in line) and is_precise:
            vcf_new.write("%s"%line)
            continue

        a = line.split()
        chr_name = a[0]
        if chr_name != chr_name_old:
            print("Processing chromosome",chr_name)
            chr_change = True
            chr_name_old = chr_name
            del fa
            fa = util.readGenome(genome,chr_name_old)

        start = int(a[1])
        tmp = a[7].split(';')
        tmp2 = tmp[0].split('=')
        end = int(tmp2[1])
        sv = a[4]
        if end>=start:
            if "DEL" in sv:
                for i in range(len(a)):
                    if i == 3:
                        vcf_new.write("%s\t"%fa[start])
                    elif i == 4:
                        vcf_new.write("%s\t"%fa[start:end])
                    else:
                        vcf_new.write("%s\t"%a[i])
                vcf_new.write("\n")
            elif "INV" in sv:
                for i in range(len(a)):
                    if i == 3:
                        vcf_new.write("%s\t"%fa[start:end])
                    elif i == 4:
                        vcf_new.write("%s\t"%reverseComplement(fa[start:end]))
                    else:
                        vcf_new.write("%s\t"%a[i])
                vcf_new.write("\n")

            elif "TANDEM" in sv:
                for i in range(len(a)):
                    if i == 3:
                        vcf_new.write("%s\t"%fa[start])
                    elif i == 4:
                        vcf_new.write("%s\t"%fa[start:end])
                    else:
                        vcf_new.write("%s\t"%a[i])
                vcf_new.write("\n")
            elif "ISP" in sv and "ISINV" in a[7]:
                for i in range(len(a)):
                    if i == 3:
                        vcf_new.write("%s\t"%fa[start])
                    elif i == 4:
                        vcf_new.write("%s\t"%reverseComplement(fa[start:end]))
                    else:
                        vcf_new.write("%s\t"%a[i])
                vcf_new.write("\n")
            elif "ISP" in sv and "ISINV" not in a[7]:
                for i in range(len(a)):
                    if i == 3:
                        vcf_new.write("%s\t"%fa[start])
                    elif i == 4:
                        vcf_new.write("%s\t"%fa[start:end])
                    else:
                        vcf_new.write("%s\t"%a[i])
                vcf_new.write("\n")
            elif "INS:ME" in sv:
                #find the mei first
                tmp = sv.split(':')
                mei_annot=tmp[2][:-1]
                mei_seq=''
                if meis[mei_annot] == "":
                    mei_start,mei_end = util.findMEISeq(mei_annot_file, mei_annot, chr_name)
                    if mei_start == -1:
                        vcf_new.write("%s"%line)
                        continue
                    else:
                        mei_seq = fa[mei_start:mei_end]
                        meis[mei_annot] = mei_seq
                else:
                    mei_seq = meis[mei_annot]

                for i in range(len(a)):
                    if i == 3:
                        vcf_new.write("%s\t"%fa[mei_start])
                    elif i == 4:
                        vcf_new.write("%s\t"%mei_seq)
                    else:
                        vcf_new.write("%s\t"%a[i])
                vcf_new.write("\n")

            else:
                vcf_new.write("%s"%line)
