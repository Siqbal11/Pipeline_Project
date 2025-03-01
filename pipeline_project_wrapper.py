#Step 1: import necessary libraries that will be used to run the wrapper
import os
import subprocess
from Bio import SeqIO

Project_Directory = f"PipeProject_SamiahIqbal" 
os.system(f"mkdir {Project_Directory}") #creates a project directory
os.chdir(Project_Directory) #changes to the project directory

logfile = "PipleineProject.log" #where we are writing our output to

os.system("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030") #download the raw data from NCBI
os.system("fasterq-dump SRR5660030") #converts the SRA file to fastq-files

os.system ("wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR005/660/SRR5660033.sralite.1") #download the raw data from NCBI
os.system("fasterq-dump SRR5660033") #converts the SRA file to fastq files


#Step 2:
os.system("bowtie2-build NC_006273.2.fasta HCMV_index") #builds the HCMV_index that will be used for mapping
os.system("bowtie2 --quiet -x HCMV_index -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S mapped_reads_1.sam") #uses the paired-end reads to create a sam file
os.system("samtools fastq mapped_reads_1.sam > filtered_mapped_reads_1.fastq") #sam file is converted to a fastq file that contains the mapped reads for the first SRR

os.system("bowtie2 --quiet -x HCMV_index -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S mapped_reads_2.sam") #uses the paired-end reads to create a sam file
os.system("samtools fastq mapped_reads_2.sam > filtered_mapped_reads_2.fastq") #sam file is converted to a fastq file that contains the mapped reads for the second SRR


def count_reads(file): #count the number of reads before and after filtering
    return subprocess.getoutput(f"grep -c '^@' {file}").strip() #all reads begin with @

before_reads_1 = count_reads("SRR5660030_1.fastq") #store number of reads for before mapping
after_reads_1 = count_reads("filtered_mapped_reads_1.fastq") #store number of reads after mapping

before_reads_2 = count_reads("SRR5660033_1.fastq")  #store number of reads for before mapping
after_reads_2 = count_reads("filtered_mapped_reads_2.fastq") #store number of reads after mapping

with open(logfile, "a") as log: #writes the output to a log file
    log.write(f"Donor1_2dpi had {before_reads_1} before Bowtie2 filtering and {after_reads_1} read pairs after.")
    log.write(f"Donor1_2dpi had {before_reads_2} before Bowtie2 filtering and {after_reads_2} read pairs after.")     

#Step 3:
os.system("nohup spades.py -k 99 -t 2 --only-assembler -1 filtered_mapped_reads_1.fastq -2 filtered_mapped_reads_2.fastq -o spades_assembly/ &") #uses 99 as the kmer and runs the spades assembly to get the contigs.fasta file

#Step 4:
def analyze_contigs(fasta_file):
    count_contigs = 0 #initialize the count
    total_bp = 0 #initialize the count
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contig_length = len(record.seq) #look a the length
            if contig_length > 1000: #if length is greater than 1000 bo, then count contigs
                count_contigs += 1 #increment count_contigs if the above statement
                total_bp += contig_length #calculates the total number of bp in the contigs that are longer than 1000 bp

    return count_contigs, total_bp

fasta_file = "spades_assembly/contigs.fasta"

with open(logfile, "a") as log: #writes output to a log file
    log.write(f"There are {count_contigs} > 1000 bp in the assembly")
    log.write(f"There are {total_bp} bp in the assembly.")

#Step 5: #function is used to find the longest contig
def longest_contig(fasta_file):
    max_length = 0 #initialize it to be 0
    with open(contigs_assembly, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            contig_length = len(record.seq) #calculates the length of a contig
            if contig_length > max_length: #if a longer contig is found
                longest_contig = record #update it to the record
                max_length = contig_length #update to the new longest contig length
    return longest_contig #return the longest contig

longest_contig_stored = longest_contig(fasta_file) #store the longest contig
with open("longest_contig.fasta", "w") as handle:
    SeqIO.write(longest_contig_stored, handle, "fasta")

os.system("makeblastdb -in betaherpesvirinae.fasta -out betaherpes_db -title Betaherpesvirinae -dbtype nucl") #make a database
os.system("blastn -query longest_contig.fasta -db betaherpes_db -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' -out blast_results.txt") #run blast to query the nr nucleotide database

with open("blast_results.txt", "r") as blast_results, open(logfile, "a") as log: #view the file that contains the blast results
    log.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n") #this writes the header in the log file
    for i, line in enumerate(blast_results): #go through each line
        if i >= 10: #10 results
            break
        log.write(line) #write each line to the log file
