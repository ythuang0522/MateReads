# Introduction
MateReads is an experimental and extended module from StriDe assembler: https://github.com/ythuang0522/StriDe. It converts mate-pair reads from long insert library into super long reads, in order to provide additional scaffolding material (long reads) other than mate-pair reads. Theoretically, long reads can help scaffold highly-fragmented contigs where mate-pair reads over-span distant contigs and lose the adjacent linkage between small contigs in proximity.

# Compilation 
MateReads depends on google sparsehash and zlib. Please install these two libraries before compilation.

      1. ./autogen.sh 
      2. ./configure
      3. make -j 8

After compilation, there will be an executable file named "stride" in the StriDe folder.

# Execution of MateReads
The MateReads require input of paired-end reads (from short insert) and mate-pair reads (from long insert). The paired-end reads are used to construct an FM-index, which will be used to convert mate-pair reads into long reads. Suppose pe.fa represents paired-end reads and mate.fa is mate-pair reads. The following commment build the FM-index of short reads (pe.fa).

      ./stride index -t 20 pe.fa

The following command converts mate-pair reads (mate.fa) into long reads (mate.merge.fa).

      ./stride matepair -M 90 -t 30 -I 3500 -c 437 -p pe mate.fa

Note that mate.fa (fasta) should be outward and interleaved in the file.

Major Arguments: 

-M : the maximum size of FM-index extension, which should be below paired-end read length (recommend PE read length * 0.9). 

-I : the mamimum insert size of long insert library.

-c : the minimum insert size of long long insert library (default: 375 or set 0.125 * maximum insert size)

-p : the prefix of FM-index of paired-end reads.

Advanced Arguments:

-x : the threshold of k-mer in each k-mer extension

-t : the threads to be used. 

-k : the kmer size for trimming kmer.

-m : the refinement size for the k-mer when resizing 

-L : the mamimum leave to be used.



----------------------------------------------------------------------------------------------

Output:

The program produce 4 files.

1. mate.merge.fa : Long read sequences are stored in this file.

2. mate.shortIS.fa: Possibly-contaminated short-insert reads.

3. mate.trimmed.fa: Mate pair reads failed to be converted into long sequences.

4. mate.polluted.fa:  One end of mate pair reads is possibly contamination (i.e., no sufficient k-mer frequency).
