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
