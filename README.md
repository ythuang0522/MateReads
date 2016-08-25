# Execution of Matepair Extension
matepair extension is modificated from FM-walk in StriDe.
StriDe: https://github.com/ythuang0522/StriDe

Usage:

stride fmwalk -m 30 -M 90 -a trimMatePair -t 30 -L 64 -I 3500 -x 3 -k 31 -p ecor inputfile.fa

inputfile.fa : Input file(fasta), matepair should be outward and interleaved in the file.

-a : the FM-walk algorithm. trimMatePair is for matepair extension

-m : the refinement size for the k-mer when resizing 

-M : the maximum size of extension k-mer. When k-mer exceeds maximum size, the k-mer resize into minimum size.

-L : the mamimum leave to be used.

-I : the mamimum length to be extened.

-x : the threshold of k-mer in each k-mer extension

-t : the threads to be used. 

-p : the prefix of FM-index file.

-k : the kmer size for trimming kmer.

.
