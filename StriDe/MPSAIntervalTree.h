//----------------------------------------------
// Copyright 2014 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------

//
// Re-written from Jared Simpson's StringThreaderNode and StringThreader class
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef MPSAIntervalTree_H
#define MPSAIntervalTree_H

#include <list>
#include "BWT.h"
#include "HashMap.h"
#include "BWTAlgorithms.h"


// Typedefs
class MPSAIntervalNode;
typedef std::list<MPSAIntervalNode*> STNodePtrList;

// Object to hold the result of the threading process
struct MPSAIntervalNodeResult
{
    std::string thread;
	size_t SAICoverage;
};
typedef std::vector<MPSAIntervalNodeResult> MPSAIntervalNodeResultVector;

// A node in the threading tree
class MPSAIntervalNode
{
    public:

        //
        // Functions
        //
        MPSAIntervalNode(const std::string* pQuery, MPSAIntervalNode* parent);
        ~MPSAIntervalNode();

        // Add a child node to this node with the given label
        // Returns a pointer to the created node
        MPSAIntervalNode* createChild(const std::string& label);

        // Extend the label of this node by l
        void extend(const std::string& ext);

        // Return a suffix of length l of the string represented by this node
        std::string getSuffix(size_t l) const;

        // Return the complete sequence of the string represented by the branch
        std::string getFullString() const;


        // Initialize or update the alignment data
        //void computeInitialAlignment(const std::string& initialLabel, int queryAlignmentEnd, int bandwidth);
        void computeInitial(const std::string& initialLabel);

        // Recursive function to print all the strings represented
        // by this node and all its children.
        void printAllStrings(const std::string& parent) const;


        size_t getKmerCount(){return m_kmerCount;};
        void addKmerCount(size_t kmercount){ m_kmerCount+=kmercount; };
        // size_t getKmerSize(){return m_kmersize;};
        // void setKmerSize(size_t kmersize){ m_kmersize=kmersize; };

        BWTInterval fwdInterval;
        BWTInterval rvcInterval;

    private:

        //
        // Data
        //

        // The extension string from the parent
        std::string m_label;
        size_t m_kmerCount;
        // size_t m_kmersize;


        // The query string being threaded through the graph
        const std::string* m_pQuery;

        // The parent node, can be NULL
        MPSAIntervalNode* m_pParent;
        STNodePtrList m_children;


};

class MPSAIntervalTree
{
    public:
        MPSAIntervalTree(const std::string* pQuery,
                       size_t minOverlap,
					   size_t maxOverlap,
                       size_t MaxLength,
                       size_t MaxLeaves,
					   BWTIndexSet indices,
                       std::string secondread,
                       //const size_t repeat_threshold,
                       size_t SA_threshold=3,//chaohung103 20151216
                       bool KmerMode=false);
        MPSAIntervalTree(const std::string* pQuery,
                       size_t minOverlap,
                       size_t maxOverlap,
                       size_t MaxLength,
                       size_t MaxLeaves,
					   BWTIndexSet indices,
                       const size_t repeatThreshold,
                       std::string secondread,
                       //const size_t repeat_threshold,
                       size_t SA_threshold=3,//chaohung103 20151216
                       bool KmerMode=false);
        ~MPSAIntervalTree();

        //return the merged string
        //bool mergeTwoReads(StringVector & mergeReads);
        int mergeTwoReads(std::string &mergedseq);
        int mergeRepeat(std::string &mergedseq); //repeat case - merging round 2 
		size_t getKmerCoverage(){return m_maxKmerCoverage;};
		size_t getMaxUsedLeaves(){return m_maxUsedLeaves;};
        size_t getCurrentExtendAllFreq(){return m_currentExtendAllFreq;};
        void addCurrentExtendAllFreq(size_t kmerfreq){ m_currentExtendAllFreq+=kmerfreq; };
        void initCurrentExtendAllFreq(void){ m_currentExtendAllFreq= 0; };
		bool isBubbleCollapsed(){return m_isBubbleCollapsed;}
        // Print all the strings represented by the tree
        void printAll();

    private:

        //
        // Functions
        //
        void extendLeaves();
        void attempToExtend(STNodePtrList &newLeaves);
        void filterLowCoverageInterval(STNodePtrList &newLeaves);
        void refineSAInterval(size_t newKmerSize);
        std::vector<std::pair<std::string, BWTIntervalPair> > getFMIndexExtensions(MPSAIntervalNode* pNode);

        // Check if the leaves can be extended no further
        bool isTerminated(MPSAIntervalNodeResultVector& results);
        bool isTwoReadsOverlap(std::string & mergedseq);
        size_t calculateKmerCoverage (const std::string & seq , size_t kmerLength , const BWT* pBWT);
		bool replaceLowFreqKmer (std::string & seq , size_t kmerLength);
        void testingOneEND(MPSAIntervalNodeResultVector& results);
        void removeLeavesByRepeatKmer();

        //
        // Data
        //
        const std::string* m_pQuery;
        size_t m_minOverlap;
		size_t m_maxOverlap;
        size_t m_MaxLength;
        size_t m_MaxLeaves;
        
		BWTIndexSet m_indices;
        const size_t m_repeatThreshold;
        std::string m_secondread;
        size_t m_min_SA_threshold;
        bool m_kmerMode;

        MPSAIntervalNode* m_pRootNode;
        STNodePtrList m_leaves;
        
        DenseHashMap<std::string, size_t, StringHasher> m_KmerIndexMap;
        size_t m_currentLength;
		size_t m_currentKmerSize;
		size_t m_maxKmerCoverage;
		size_t m_maxUsedLeaves;
        bool m_isBubbleCollapsed;
        size_t m_currentExtendAllFreq;
        

        
        BWTInterval m_fwdTerminatedInterval;   //in rBWT
        BWTInterval m_rvcTerminatedInterval;   //in BWT
};

#endif
