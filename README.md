# VCFresolve

VCFresolve aims to sequence resolve the VCF files...

# Downloading VCFresolve
	git clone https://github.com/asylvz/VCFresolve.git

# Quick Start (To resolve the sequences of the SVs)
    python vcf-resolve --resolve -g Reference_Genome -v VCF -a Repeat_Annotation -o Output
# Optionally you can sequence resolve only the precise calls, i.e., the calls that don't have IMPRECISE keyword
    python vcf-resolve --resolve -g Reference_Genome -v VCF -a Repeat_Annotation -o Output --precise
    
# Running Mendelian Filter for trios
    python vcf-resolve --mendelian -v VCF -o Output


# Downloading the repeat annotation file

Using UCSC sequence and annotations [download page](http://hgdownload.cse.ucsc.edu/downloads.html);

For repeat annotations, navigate to the related genome's "Full data set" page and download the RepeatMasker .out files (i.e., http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromOut.tar.gz is for human genome GRCh37/hg19)

Extract the files into a folder and merge them into a single .out file:
	
	cat * >reps.out
	
*Alternatively, if the files are inside directories for each chromosome, then use:

	cat */* >reps.out
	

# Downloading the reference genome file

Using UCSC sequence and annotations [download page](http://hgdownload.cse.ucsc.edu/downloads.html);

For reference genome, navigate to the related genome's "Full data set" page and download the FASTA files (i.e., http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz is for human genome GRCh37/hg19)

Extract the files into a folder and merge them into a single .fasta file:
	
	cat * >ref.fasta
