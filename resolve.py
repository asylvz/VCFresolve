import os
import utilities as util

def resolve(genome, vcf_file, vcf_new, mei_annot_file, is_precise, vcf_folder, verbose):

    vcfs=[]
    if not vcf_folder is None and not vcf_folder[-1] == '/':
        vcf_folder +="/"

    print()
    if is_precise:
        print("Sequence resolving the VCF file for only the precise calls")
    else:
        print("Sequence resolving the VCF file for all the calls")

    file_count = 0
    if not vcf_folder is None:
        for filename in os.listdir(vcf_folder):
            if filename.endswith(".vcf"):
                vcfs.append(vcf_folder+filename)
                file_count +=1
    else:
        vcfs.append(vcf_file)
        file_count +=1

    for i in range(len(vcfs)):
        if not vcf_folder is None:
            tmp = vcfs[i].split('.vcf')
            vcf_new = vcf_new = open(tmp[0]+"_resolved.vcf","w")

        vcf, comments = util.readVcf(vcfs[i])
        #chr_name_old = "1"
        #fa = util.readGenome(genome,chr_name_old)
        beginning = True
        for line in comments:
            vcf_new.write("%s"%line)

        meis = dict()

        for line in vcf:
            a = line.split()
            if "INS:ME" in a[4]:
                tmp = a[4].split(':')
                mei_annot=tmp[2][:-1]
                meis.update({mei_annot:""})

        print("Running for",vcfs[i],"(",i+1,"of",file_count,")")

        fa = ''
        for line in vcf:
            if ("IMPRECISE" in line) and is_precise:
                vcf_new.write("%s"%line)
                continue


            a = line.split()
            chr_name = a[0]
            if beginning:
                chr_name_old = chr_name
            if beginning or chr_name != chr_name_old:
                if verbose:
                    print("Processing chromosome", chr_name)
                beginning = False
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
                            vcf_new.write("%s\t"%fa[start:end])
                        elif i == 4:
                            vcf_new.write("%s\t"%fa[start])
                        else:
                            vcf_new.write("%s\t"%a[i])
                    vcf_new.write("\n")
                elif "INV" in sv:
                    for i in range(len(a)):
                        if i == 3:
                            vcf_new.write("%s\t"%fa[start:end])
                        elif i == 4:
                            vcf_new.write("%s\t"%util.reverseComplement(fa[start:end]))
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
        vcf_new.close()
