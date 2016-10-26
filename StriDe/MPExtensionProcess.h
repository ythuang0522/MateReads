//----------------------------------------------
// Copyright 2014 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// MPExtensionProcess - Wrapper to perform error correction
// for a sequence work item
//
#ifndef MPExtensionProcess_H
#define MPExtensionProcess_H

#include "HashMap.h"
#include "Util.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "Metrics.h"
#include "BWTIndexSet.h"
#include "SampledSuffixArray.h"
#include "BWTAlgorithms.h"
#include "BitVector.h"
#include "KmerDistribution.h"

enum MPExtensionAlgorithm
{

    //trimmatepair by chaohung 20151017
    MPE_trimMatePair
};


enum NextKmerDir
{
	NK_START,
	NK_END
};

// Parameter object for the error corrector
struct MPExtensionParameters
{
    BWTIndexSet indices;

    int numKmerRounds;
    int kmerLength;

    // output options
    bool printOverlaps;

	int maxLeaves;
	int maxInsertSize;
	int minInsertSize;
    int minOverlap;
	int maxOverlap;
	size_t FreqThreshold;
	KmerDistribution	 kd;

};


struct KmerContext
{
	public:

	//empty
	KmerContext()
	{
		kmerLength=0;
		readLength=0;
		numKmer=0;
	}

	//originalSeq
	KmerContext(std::string seq,size_t kl, BWTIndexSet & index)
	{
		if( seq.length() >= kl)
		{
			readSeq = seq ;
			readLength = readSeq.length();
			kmerLength = kl ;
			numKmer = readLength-kmerLength+1 ;
			kmers.resize(numKmer);
			kmerFreqs_same.resize(numKmer);
			kmerFreqs_revc.resize(numKmer);

			for (size_t i = 0 ; i < numKmer  ; i++)
			{
				kmers[i] = readSeq.substr (i,kmerLength);
				kmerFreqs_same[i] = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmers[i], index) ;
				kmerFreqs_revc[i] = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmers[i]), index) ;
			}
		}
		else
		{
			kmerLength=0;
			readLength=0;
			numKmer=0;
		}
	}

	//subSeq
	KmerContext( KmerContext origin , int head , int tail)
	{
		assert (head>=0 && tail>=0 && head<=tail) ;
		kmerLength = origin.kmerLength;
		readSeq= origin.readSeq.substr(head,tail-head+kmerLength);
		readLength = readSeq.length();
		kmers.assign (origin.kmers.begin()+head,origin.kmers.begin()+tail+1);
		kmerFreqs_same.assign (origin.kmerFreqs_same.begin()+head , origin.kmerFreqs_same.begin()+tail+1);
		kmerFreqs_revc.assign (origin.kmerFreqs_revc.begin()+head , origin.kmerFreqs_revc.begin()+tail+1);

		assert (kmers.size() == kmerFreqs_same.size());
		assert (kmers.size() == kmerFreqs_revc.size());
		assert (readLength-kmerLength+1 == kmers.size());
		numKmer = kmers.size();
	}


	std::string readSeq;
	size_t kmerLength;

	size_t readLength;
	size_t numKmer ;
	std::vector<std::string> kmers;
	std::vector<size_t> kmerFreqs_same;
	std::vector<size_t> kmerFreqs_revc;

	bool empty(){ return readSeq.empty() ;}

};

class MPExtensionResult
{
    public:
        MPExtensionResult()
		: kmerize(false),kmerize2(false),merge(false),merge2(false),trim(false),polluted(false), mergeMP(false), tryMerge(false),longmergeCount(false),shortmergeCount(false),whyLongFail(0), head1(false),head2(false){}

        DNAString correctSequence;
		DNAString correctSequence2;
        bool kmerize;
		bool kmerize2;
		bool merge;
		bool merge2;
        size_t kmerLength;
		std::vector<DNAString> kmerizedReads ;
		std::vector<DNAString> kmerizedReads2 ;
        //trimmatepair by chaohung 20151017
		DNAString pollutedSequence;
		DNAString pollutedSequence2;
        bool trim;
		bool polluted;
		bool mergeMP;
		bool tryMerge;
		bool longmergeCount;
		bool shortmergeCount;
		int whyLongFail;
        bool head1;
        bool head2;

};

//
class MPExtensionProcess
{
    public:
        MPExtensionProcess(const MPExtensionParameters params);
        ~MPExtensionProcess();

        MPExtensionResult process(const SequenceWorkItem& item);
        MPExtensionResult correct(const SequenceWorkItem& item);

		MPExtensionResult kmerTrimCorrection(const SequenceWorkItem& workItem);
		MPExtensionResult kmerizeLowKmerReadCorrection(const SequenceWorkItem& workItem);

		/***************************************************************************/

		MPExtensionResult process(const SequenceWorkItemPair& workItemPair)
		{
			// return mergePairEndCorrection(workItemPair);
				//trimmatepair by chaohung 20151017

			return MatepairWorkflow(workItemPair);
			MPExtensionResult result;
			return result;
		}
		
        MPExtensionResult MatepairWorkflow(const SequenceWorkItemPair& workItemPair);//trimmatepair by chaohung 20151017

    private:		
		//check necessary conditions for FM-index walk
		
		std::string getReliableInterval(std::string& seq, KmerContext& kc);

		size_t numNextKmer(std::string kmer , NextKmerDir dir ,BWTIndexSet & index, size_t threshold);
		bool isSimple (std::string Lkmer, std::string Rkmer, BWTIndexSet & index, size_t threshold) ;

		bool existStrongLink (std::string Lkmer,std::string Rkmer,BWTIndexSet & index,size_t threshold) ;
		bool existNextStrongKmer(std::string kmer , NextKmerDir dir ,BWTIndexSet & index,size_t threshold) ;
		bool isIntervalExistStrongKmer (std::pair<size_t,size_t> interval,std::vector<size_t> & countQualified);
		bool isPathReliable(std::pair<size_t,size_t> intervalX, std::pair<size_t,size_t> intervalY,std::vector<size_t> & countQualified);
		bool isIntervalMerge (std::vector< std::pair<size_t,size_t> > & intervals , std::vector<size_t> & countQualified );

		//trim dead-end by de Bruijn graph using FM-index
        std::string trimRead ( std::string readSeq ,size_t kmerLength ,size_t threshold ,BWTIndexSet & index);
        //trimmatepair by chaohung 20151017
		//std::string trimRead ( std::string readSeq ,size_t kmerLength ,size_t threshold ,BWTIndexSet & index,bool *H, int *cs1, int *a1, const int freqThrshold);
        std::string trimRead( std::string readSeq ,size_t kmerLength , size_t /*threshold*/ ,BWTIndexSet & index,bool *H, const size_t freqThrshold);
		int splitRead (KmerContext& seq, std::vector<std::string> & kmerReads ,size_t threshold, BWTIndexSet & index);

		// bool hasPESupport (std::string r1,std::string r2
	                     // , BWTIndexSet & index , ReadInfoTable*  pRIT
						 // , size_t firstK , size_t secondK);


		int getMainSeed (KmerContext seq, std::vector<KmerContext> & kmerReads ,size_t threshold,BWTIndexSet & index);
		//split read to kmers
		std::vector<size_t> splitRead( KmerContext seq ,size_t threshold ,BWTIndexSet & index ,size_t singleThreshld =0 );

		bool  isLowComplexity (std::string seq , float & GCratio);
		size_t maxCon (std::string s);

        MPExtensionParameters m_params;




};

// Write the results from the overlap step to an ASQG file
class MPExtensionPostProcess
{
    public:
        MPExtensionPostProcess(std::ostream* pCorrectedWriter,
                                std::ostream* pDiscardWriter,
                                const MPExtensionParameters params);
		//trimmatepair by chaohung 20151030
		MPExtensionPostProcess(std::ostream* pCorrectedWriter,
                                std::ostream* pDiscardWriter,
								std::ostream* pLongReadWriter,
								std::ostream* pShortReadWriter,
                                const MPExtensionParameters params);
        ~MPExtensionPostProcess();

        void process(const SequenceWorkItem& item, const MPExtensionResult& result);
		void process(const SequenceWorkItemPair& itemPair, const MPExtensionResult& result);

    private:

        std::ostream* m_pCorrectedWriter;
        std::ostream* m_pDiscardWriter;
		std::ostream* m_pLongReadWriter;
		std::ostream* m_pShortReadWriter;
        std::ostream* m_ptmpWriter;
		MPExtensionParameters m_params;
        // DenseHashSet<std::string,StringHasher> *m_pCachedRead;

		size_t m_kmerizePassed ;
		size_t m_mergePassed ;
        size_t m_qcFail;
        //trimmatepair by chaohung 20151017
        size_t m_trimPassed ;
        size_t m_r1_pass;
        size_t m_r2_pass;
		size_t m_totalMerge;
		size_t m_longSucc;
		size_t m_HighError;
		size_t m_exceedSearchDepth;
		size_t m_Repeat;
		size_t m_Case4;


};

#endif
