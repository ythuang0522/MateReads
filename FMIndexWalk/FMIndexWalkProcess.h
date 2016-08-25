//----------------------------------------------
// Copyright 2014 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// FMIndexWalkProcess - Wrapper to perform error correction
// for a sequence work item
//
#ifndef FMIndexWalkProcess_H
#define FMIndexWalkProcess_H

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

enum FMIndexWalkAlgorithm
{
	FMW_KMERIZE,
	FMW_MERGE,
	FMW_HYBRID,
    //trimmatepair by chaohung 20151017
    FMW_trimMatePair
};


enum NextKmerDir
{
	NK_START,
	NK_END
};

// Parameter object for the error corrector
struct FMIndexWalkParameters
{
    FMIndexWalkAlgorithm algorithm;
    BWTIndexSet indices;

    int numKmerRounds;
    int kmerLength;

    // output options
    bool printOverlaps;

	int maxLeaves;
	int maxInsertSize;
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

class FMIndexWalkResult
{
    public:
        FMIndexWalkResult()
		: kmerize(false),kmerize2(false),merge(false),merge2(false),trim(false),polluted(false), mergeMP(false), tryMerge(false),longmergeCount(false),shortmergeCount(false),whyLongFail(0),whyShortFail(0), head1(false),head2(false){}//,c1(0),p1(0){}

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
		//int cs1[100];
		//int ps1[100];
        bool trim;//By Chaohung_2015_01_22 used to control output
		bool polluted;//By Chaohung_2015_03_10 used to control output
		bool mergeMP;//By Chaohung_2015_10_30 used to control output
		bool tryMerge;
		bool longmergeCount;
		bool shortmergeCount;
		int whyLongFail;
		int whyShortFail;
        bool head1;
        bool head2;
		//int c1;
		//int p1;
		

};

//
class FMIndexWalkProcess
{
    public:
        FMIndexWalkProcess(const FMIndexWalkParameters params);
        ~FMIndexWalkProcess();

        FMIndexWalkResult process(const SequenceWorkItem& item);
        FMIndexWalkResult correct(const SequenceWorkItem& item);

		FMIndexWalkResult kmerTrimCorrection(const SequenceWorkItem& workItem);
		FMIndexWalkResult kmerizeLowKmerReadCorrection(const SequenceWorkItem& workItem);

		/***************************************************************************/

		FMIndexWalkResult process(const SequenceWorkItemPair& workItemPair)
		{
			// return mergePairEndCorrection(workItemPair);
			switch(m_params.algorithm)
			{
				case FMW_HYBRID:
					{
						return MergeAndKmerize(workItemPair);
						break;
					}
				case FMW_MERGE:
				{
					return MergePairedReads(workItemPair);
					break;
				}
                //trimmatepair by chaohung 20151017
                case FMW_trimMatePair:
				{
					return MatepairWorkflow(workItemPair);
					break;
				}
				default:
				{
						assert(false);
				}
			}
			FMIndexWalkResult result;
			return result;
		}
		
		FMIndexWalkResult MergeAndKmerize(const SequenceWorkItemPair& workItemPair);
		FMIndexWalkResult MergePairedReads(const SequenceWorkItemPair& workItemPair);
        FMIndexWalkResult MatepairWorkflow(const SequenceWorkItemPair& workItemPair);//trimmatepair by chaohung 20151017

    private:		
		//check necessary conditions for FM-index walk
		bool isSuitableForFMWalk(std::string& seqFirst, std::string& seqSecond);
		
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

        FMIndexWalkParameters m_params;




};

// Write the results from the overlap step to an ASQG file
class FMIndexWalkPostProcess
{
    public:
        FMIndexWalkPostProcess(std::ostream* pCorrectedWriter,
                                std::ostream* pDiscardWriter,
                                const FMIndexWalkParameters params);
		//trimmatepair by chaohung 20151030
		FMIndexWalkPostProcess(std::ostream* pCorrectedWriter,
                                std::ostream* pDiscardWriter,
								std::ostream* pLongReadWriter,
								std::ostream* pShortReadWriter,
                                const FMIndexWalkParameters params);
        ~FMIndexWalkPostProcess();

        void process(const SequenceWorkItem& item, const FMIndexWalkResult& result);
		void process(const SequenceWorkItemPair& itemPair, const FMIndexWalkResult& result);

    private:

        std::ostream* m_pCorrectedWriter;
        std::ostream* m_pDiscardWriter;
		std::ostream* m_pLongReadWriter;
		std::ostream* m_pShortReadWriter;
        std::ostream* m_ptmpWriter;
		FMIndexWalkParameters m_params;
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
		size_t m_shortSucc;
		size_t m_HighError;
		size_t m_exceedSearchDepth;
		size_t m_Repeat;
		size_t m_Case4;
		size_t m_HighError2;
		size_t m_exceedSearchDepth2;
		size_t m_Repeat2;
		size_t m_Case42;
        //size_t s_numHighError;
		//size_t s_numExceedSearchDepth;
		//size_t s_numRepeat;
		//size_t s_numCase4;
        //size_t s_numWalkTotal;
        //size_t s_numWalkSucc;
        //size_t s_numWalkFail;
};

#endif
