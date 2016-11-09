# Introduction
MateReads is an experimental and extended module from StriDe assembler: https://github.com/ythuang0522/StriDe. It converts mate-pair reads from long insert library into super long reads, in order to provide additional scaffolding material (long reads) other than mate-pair reads. Theoretically, long reads can help scaffold highly-fragmented contigs where mate-pair reads over-span distant contigs and lose the adjacent linkage between small contigs in proximity.

# To compile in your environment, type 

      1. ./autogen.sh 
      2. ./configure
      3. make

and there will be an executable file named stride in the StriDe folder.

# Execution of MateReads
The MateReads require input of paired-end reads (short insert) and mate-pair reads (long insert). The paired-end reads are used to construct an FM-index, which will be used to convert mate-pair reads into long reads. Suppose pe.fa represents paired-end reads and mate.fa is mate-pair reads. The following commment build the FM-index of short reads.

      ./stride index -t 20 pe.fa

The following command converts mate-pair reads into long reads.

      ./stride matepair -m 30 -M 90 -t 30 -L 64 -I 3500 -c 437 -x 3 -k 31 -p pe mate.fa

----------------------------------------------------------------------------------------------
mate.fa : Input file(fasta), matepair should be outward and interleaved in the file.
----------------------------------------------------------------------------------------------

Argument : 

-m : the refinement size for the k-mer when resizing 

-M : the maximum size of extension k-mer. When k-mer exceeds maximum size, the k-mer resize into minimum size.

-L : the mamimum leave to be used.

-I : the mamimum insert size to be extened.

-c : the minimum insert size of long sequence is to be reserved (default: 375),usually set 0.125 * maximum insert size

-x : the threshold of k-mer in each k-mer extension

-t : the threads to be used. 

-p : the prefix of FM-index of paired-end reads.

-k : the kmer size for trimming kmer.

----------------------------------------------------------------------------------------------

Output:

The program produce 4 files.

1. inputfile.merge.fa : Long read sequences are reserved in this file.

2. inputfile.shortIS.fa: Long read sequences, which lengths are shorter than the minimum insert size (-c option), are reserved in the file.

3. inputfile.trimmed.fa: Mate pair reads, which can't be converted into long sequences, are reserved in the file.

4. inputfile.polluted.fa:  one of the end of mate pair reads can not be found sufficient k-mer in FM-index.
