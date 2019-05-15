import sys

if __name__ == "__main__":
    import argparse
    import mendelian as mf
    import resolve as rslv

    parser = argparse.ArgumentParser(description="VCFresolve")
    #parser.add_argument("echo")
    parser.add_argument("-g","--genome", dest="genome_file",
                        help="The reference genome file in fasta format.")
    parser.add_argument("-v","--vcf", dest="vcf_file", required = True,
                        help="The VCF file to resolve the sequences.")
    parser.add_argument("-a","--mei", dest="annot_file",
                        help="The Repeat Annotation files.")
    parser.add_argument("-o","--output", dest="output_file",
                        help="The output VCF file.")
    parser.add_argument("--precise", action='store_true',
                        help="To resolve the sequences of only the precise calls. The default is to resolve all the calls.")
    parser.add_argument("--resolve", action='store_true',
                        help="Used to find the sequences of the SVs.")
    parser.add_argument("--mendelian", action='store_true',
                        help="This option allows you to use mendelian filter for trios. With this option, only the VCF file (-v, --vcf) and output file (-o, --output) is needed as input.")
    args = parser.parse_args()


    #print("The arguments given are: ")
    #print("Genome is",args.genome_file)
    #print("VCF File is", args.vcf_file)
    #print("MEI Annotation file is", args.annot_file)
    #print("Output File is", args.output_file)7

    if args.output_file is None:
        output = args.vcf_file+"_v2"
        vcf_new = open(output,"w")
    else:
        output = args.vcf_file
        vcf_new = open(output,"w")

    if args.resolve:
        if args.annot_file is None or args.genome_file is None:
            sys.exit('Please run with correct genome and annotation file')
        else:
            mei_annot_file = args.annot_file
            genome = args.genome_file
            rslv.resolve(genome, args.vcf_file, vcf_new, mei_annot_file, args.precise)

    elif args.mendelian:
        print("Running mendelian filter (Make sure that your VCF has the required columns for the trio and son must be at the last column)")
        mf.flt(args.vcf_file, vcf_new)
