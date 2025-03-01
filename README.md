Step 1: 
I created a directory based on the instructions, and then I changed into that directory so that the code is in that directory. 
I set up a log file, where my output will go. 
I used the links from the handout to get the SRR numbers. I first clicked on the link and then I clicked the hyperlink under the SRR number which led me to the 5 tabs, I clicked the data access tab.
I copied the first link.
I used wget and that link to download the raw data from NCBI.
I used fasterq-dump to uncompress the data. 
I did this for both SRR numbers.

Step 2:
I first built the index- HCMV_index
I used the paired-end reads to create a sam file
I then converted the sam file to a fastq file
I did this for the second SRR number
I used the count_reads function to count the before and after filtering reads for both SRR numbers
I wrote the output to the log file.

Step 3: 
I used spades with 99 as the kmer and I created a spades_assembly

Step 4: 
I used a function to analyze the contigs and the total base pairs in the assembly using the contigs.fasta file from the assembly.
If the length was greater than 1000 bp, then the count_contigs was incremented
Gives us the total base pairs based on the longest contig
I wrote that to the output 

Step 5:
I used another function to find the longest contig
I made a database
I ran blast and I used a txt file
I wrote the top 10 hits


