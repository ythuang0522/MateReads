# Execution of Matepair Extension
matepair extension is modificated from FM-walk in StriDe.
StriDe: https://github.com/ythuang0522/StriDe

Usage:

./stride matepair -m 30 -M 90 -t 30 -L 64 -I 3500 -x 3 -k 31 -p ecor inputfile.fa

inputfile.fa : Input file(fasta), matepair should be outward and interleaved in the file.
	eg.	In inputfile.fa, the sequences are ordered as follows.
		>100/1
		ACCGTCCTGTACTACT
		>100/2
		GCGATCGAACGGGACT
-m : the refinement size for the k-mer when resizing 

-M : the maximum size of extension k-mer. When k-mer exceeds maximum size, the k-mer resize into minimum size.

-L : the mamimum leave to be used.

-I : the mamimum insert size to be extened.

-c : the minimum insert size of long sequence is to be reserved (default: 375),usually set 0.125 * maximum insert size

-x : the threshold of k-mer in each k-mer extension

-t : the threads to be used. 

-p : the prefix of FM-index file.

-k : the kmer size for trimming kmer.


Output:

The program produce 4 files.
-inputfile.merge.fa : Long read sequences are reserved in this file.

-inputfile.shortIS.fa: Long read sequences, which lengths are shorter than the minimum insert size (-c option), are reserved in the file.

-inputfile.trimmed.fa: Mate pair reads, which can't be converted into long sequences, are reserved in the file.
eg. In .trimmed.fa, the sequences are ordered as follows.
>100/1 -original read
ATCGGGACGGGTTTTGGGGAGAGAGAGA
>100/1 -trimmed read
ATCGGGACGGGTT
>100/2 -original read
GGCGAAAATTTCGACGACGAAACGGTTT
>100/2 -trimmed read
GGCGAAAATTTCGACGAC

-inputfile.polluted.fa:  one of the end of mate pair reads can not be found sufficient k-mer in FM-index.
