
# coding: utf-8

# In[2]:


def readGenome(filename, chr_name):
    genome =''
    k = 0
    chr_check = False
    skip = False
    chr_fq=''
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '>':
                chr_fq = line[1:].rstrip()
            else:
                if chr_fq == chr_name:
                    genome += line.rstrip()
    return genome


# In[3]:


def readVcf(filename):
    vcf=[]
    with open(filename, 'r') as f:
        for line in f:
            if line[0] != '#':
                vcf.append(line)
    return vcf


# In[35]:


def findMEISeq(filename, mei):
    with open(filename, 'r') as f:
        k = 0
        for line in f:
            if k>2:
                a=line.split();
                if a[8] == "C":
                    continue
                mei_annot = a[9]
                if mei_annot == mei:
                    chr_name = a[4]
                    start = int(a[5])
                    end = int(a[6])
                    return start,end
            else:
                k +=1
        return NULL


# In[7]:


import sys

def findMEISeq(filename, mei, chr_name):
    with open(filename, 'r') as f:
        k = 0
        for line in f:
            if line[0] == '#':
                continue
            if k>2:
                a=line.split();
                if not len(a) > 0:
                    continue
                if a[8] == "C":
                    continue
                mei_annot = a[9]
                chr_name_annot = a[4]
                #print(mei_annot, chr_name_annot, " - ",mei, "chr"+chr_name )
                if mei_annot == mei and (chr_name_annot == "chr"+chr_name or chr_name_annot == chr_name):
                    start = int(a[5])
                    end = int(a[6])
                    return start,end
            else:
                k +=1
        return -1,-1
    
def reverseComplement(s):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t
    
def readGenome(filename, chr_name):
    genome =''
    k = 0
    chr_check = False
    skip = False
    chr_fq=''
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '>':
                chr_fq = line[1:].rstrip()
            else:
                if chr_fq == chr_name:
                    genome += line.rstrip()
    return genome

def readVcf(filename):
    vcf=[]
    comments = []
    with open(filename, 'r') as f:
        for line in f:
            if line[0] != '#':
                vcf.append(line)
            else:
                comments.append(line)
    return vcf,comments

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="VCFresolve")
    #parser.add_argument("echo")
    parser.add_argument("-g","--genome", dest="genome_file",
                        help="The reference genome file in fasta format.")
    parser.add_argument("-v","--vcf", dest="vcf_file",
                        help="The VCF file to resolve the sequences.")
    parser.add_argument("-a","--mei", dest="annot_file",
                        help="The Repeat Annotation file.")
    parser.add_argument("-o","--output", dest="output_file",
                        help="The output VCF file.")
    parser.add_argument("--precise", action='store_true',
                        help="To resolve the sequences of only the precise calls. The default is to resolve all the calls.")
    args = parser.parse_args()
    
    if args.precise:
        print("Sequence resolving the VCF file for only the precise calls")
    else:
        print("Sequence resolving the VCF file for all the calls")
    
    print("The arguments given are: ")
    print("Genome is",args.genome_file)
    print("VCF File is", args.vcf_file)
    print("MEI Annotation file is", args.annot_file)
    print("Output File is", args.output_file)
    
    mei_annot_file = args.annot_file
    vcf_new = open(args.output_file,"w")
    vcf, comments = readVcf(args.vcf_file)
    genome = args.genome_file
    
    chr_name_old = "1"
    chr_change=False
    fa = readGenome(genome,chr_name_old)
    
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
        if ("IMPRECISE" in line) and args.precise:
            vcf_new.write("%s"%line)
            continue
        
        a = line.split()
        chr_name = a[0]
        if chr_name != chr_name_old:
            print("Processing chromosome",chr_name)
            chr_change = True
            chr_name_old = chr_name
            del fa
            fa = readGenome(genome,chr_name_old)

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
                    mei_start,mei_end = findMEISeq(mei_annot_file, mei_annot, chr_name)
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
                


# In[8]:




