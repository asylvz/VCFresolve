# VCFresolve

VCFresolve aims to sequence resolve the VCF files

# Usage
    python vcf-resovlve -g Reference_Genome -v VCF -a Annotation -o Output

# Downloading the reference genome and the repeat annotations

Using UCSC sequence and annotations [download page](http://hgdownload.cse.ucsc.edu/downloads.html);

1 - For reference genome, navigate to the related genome's "Full data set" page and download the FASTA files (i.e., http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz is for human genome GRCh37/hg19)

Extract the files into a folder and merge them into a single .fasta file:
	
	cat * >ref.fasta

2 - For repeat annotations, navigate to the related genome's "Full data set" page and download the RepeatMasker .out files (i.e., http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromOut.tar.gz is for human genome GRCh37/hg19)

Extract the files into a folder and merge them into a single .out file:
	
	cat * >reps.out
	
*Alternatively, if the files are inside directories for each chromosome, then use:

	cat */* >reps.out
