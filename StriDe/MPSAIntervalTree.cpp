///----------------------------------------------
// Copyright 2014 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// MPSAIntervalTree - Iteratively construct a
// string representing a walk through an assembly graph
// matching a query sequence.
//
// The assembly graph is abstractly represented as
// an FM-index.
//
#include "MPSAIntervalTree.h"
#include "BWTAlgorithms.h"

//
// MPSAIntervalNode
//
MPSAIntervalNode::MPSAIntervalNode(const std::string* pQuery, MPSAIntervalNode* parent) : 
									   m_kmerCount(0), m_pQuery(pQuery),m_pParent(parent)
{

}

// Destructor, recurisvely delete the children of the node
MPSAIntervalNode::~MPSAIntervalNode()
{
    // Delete children
    for(STNodePtrList::iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
        delete *iter;

}

// Return a suffix of length l of the path from the root to this node
std::string MPSAIntervalNode::getSuffix(size_t l) const
{
    size_t n = m_label.size();
    if(l <= n)
    {
        return m_label.substr(n - l, l);
    }
    else
    {
        assert(m_pParent != NULL);
        return m_pParent->getSuffix(l - n) + m_label;
    }
}

// Return the full string of the path from the root to this node
std::string MPSAIntervalNode::getFullString() const
{
    if(m_pParent == NULL)
        return m_label;
    else
        return m_pParent->getFullString() + m_label;
}

// Create a new child node with the given label. Returns a pointer to the new node.
MPSAIntervalNode* MPSAIntervalNode::createChild(const std::string& label)
{
    MPSAIntervalNode* pAdded = new MPSAIntervalNode(m_pQuery, this);
    m_children.push_back(pAdded);

    //assert(!m_alignmentColumns.empty());
    //pAdded->computeExtendedAlignment(label, m_alignmentColumns.back());
    pAdded->extend(label);

    return pAdded;
}

// Extend the label of this node
void MPSAIntervalNode::extend(const std::string& ext)
{
    assert(!ext.empty());
    //assert(!m_alignmentColumns.empty());
    m_label.append(ext);
}


void MPSAIntervalNode::computeInitial(const std::string& initialLabel)
{
    m_label = initialLabel;

}


// Print the string(s) represented by this node and its children
void MPSAIntervalNode::printAllStrings(const std::string& parent) const
{
    if(m_children.empty())
    {
        std::cout << ">\n" << parent + m_label << "\n";
    }
    else
    {
        for(STNodePtrList::const_iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
            (*iter)->printAllStrings(parent + m_label);
    }
}

//
// Class: MPSAIntervalTree
MPSAIntervalTree::MPSAIntervalTree(const std::string* pQuery,
                               size_t minOverlap,
							   size_t maxOverlap,
                               size_t MaxLength,
                               size_t MaxLeaves,
                               BWTIndexSet indices,
                               std::string secondread,
                               size_t SA_threshold,
                               bool KmerMode) :
                               m_pQuery(pQuery), m_minOverlap(minOverlap), m_maxOverlap(maxOverlap), m_MaxLength(MaxLength),
                               m_MaxLeaves(MaxLeaves), m_indices(indices), m_repeatThreshold(0),
                               m_secondread(secondread), m_min_SA_threshold(SA_threshold),
                                m_kmerMode(KmerMode), m_maxKmerCoverage(0), m_maxUsedLeaves(0), m_isBubbleCollapsed(false), m_currentExtendAllFreq(0)
{
    // Create the root node containing the seed string
    m_pRootNode = new MPSAIntervalNode(pQuery, NULL);
    m_pRootNode->computeInitial(*pQuery);   //store initial str of root
    m_leaves.push_back(m_pRootNode);

    m_currentLength=pQuery->length();
	m_currentKmerSize=m_minOverlap;

	//beginning kmer is a suffix of first read
    //initialize the beginning kmer SA intervals with kmer length=m_minOverlap
    std::string beginningkmer=pQuery->substr(m_currentLength-m_minOverlap);
     // std::cout<<"beginningkmer:"<<beginningkmer<<std::endl;
     // std::cout<<"beginningkmer(reverse):"<<reverse(beginningkmer)<<std::endl;
     // std::cout<<"beginningkmer(rc):"<<reverseComplement(beginningkmer)<<std::endl;
    m_pRootNode->fwdInterval=BWTAlgorithms::findInterval( m_indices.pRBWT, reverse(beginningkmer));
    m_pRootNode->rvcInterval=BWTAlgorithms::findInterval( m_indices.pBWT, reverseComplement(beginningkmer));
    
	//ending kmer is a prefix of second read
    //initialize the ending SA intervals with kmer length=m_minOverlap
    std::string endingkmer=secondread.substr(0,m_minOverlap);
    m_fwdTerminatedInterval=BWTAlgorithms::findInterval( m_indices.pRBWT, reverse(endingkmer));
    m_rvcTerminatedInterval=BWTAlgorithms::findInterval( m_indices.pBWT, reverseComplement(endingkmer));
    
	//std::cout << m_minOverlap << ":" << beginningkmer << ":" << endingkmer << "\n";
	//getchar();
}
//
// Class: MPSAIntervalTree
MPSAIntervalTree::MPSAIntervalTree(const std::string* pQuery,
                               size_t minOverlap,
							   size_t maxOverlap,
                               size_t MaxLength,
                               size_t MaxLeaves,
                               BWTIndexSet indices,
                               const size_t repeatThreshold,
                               std::string secondread,                               
                               size_t SA_threshold,
                               bool KmerMode) :
                               m_pQuery(pQuery), m_minOverlap(minOverlap), m_maxOverlap(maxOverlap), m_MaxLength(MaxLength),
                               m_MaxLeaves(MaxLeaves), m_indices(indices), 
                               m_repeatThreshold(repeatThreshold),
                               m_secondread(secondread), m_min_SA_threshold(SA_threshold),
                                m_kmerMode(KmerMode), m_maxKmerCoverage(0), m_maxUsedLeaves(0), m_isBubbleCollapsed(false), m_currentExtendAllFreq(0)
{
    // Create the root node containing the seed string
    m_pRootNode = new MPSAIntervalNode(pQuery, NULL);
    m_pRootNode->computeInitial(*pQuery);   //store initial str of root
    m_leaves.push_back(m_pRootNode);

    m_currentLength=pQuery->length();
	m_currentKmerSize=m_minOverlap;

	//beginning kmer is a suffix of first read
    //initialize the beginning kmer SA intervals with kmer length=m_minOverlap
    std::string beginningkmer=pQuery->substr(m_currentLength-m_minOverlap);
     // std::cout<<"beginningkmer:"<<beginningkmer<<std::endl;
     // std::cout<<"beginningkmer(reverse):"<<reverse(beginningkmer)<<std::endl;
     // std::cout<<"beginningkmer(rc):"<<reverseComplement(beginningkmer)<<std::endl;
    m_pRootNode->fwdInterval=BWTAlgorithms::findInterval( m_indices.pRBWT, reverse(beginningkmer));
    m_pRootNode->rvcInterval=BWTAlgorithms::findInterval( m_indices.pBWT, reverseComplement(beginningkmer));
    
	//ending kmer is a prefix of second read
    //initialize the ending SA intervals with kmer length=m_minOverlap
    std::string endingkmer=secondread.substr(0,m_minOverlap);
    m_fwdTerminatedInterval=BWTAlgorithms::findInterval( m_indices.pRBWT, reverse(endingkmer));
    m_rvcTerminatedInterval=BWTAlgorithms::findInterval( m_indices.pBWT, reverseComplement(endingkmer));
    
	//std::cout << m_minOverlap << ":" << beginningkmer << ":" << endingkmer << "\n";
	//getchar();
}


//
MPSAIntervalTree::~MPSAIntervalTree()
{
    // Recursively destroy the tree
    delete m_pRootNode;
}

//On success return the length of merged string
int MPSAIntervalTree::mergeTwoReads(std::string &mergedseq)
{
    MPSAIntervalNodeResultVector results;
    //MPSAIntervalNodeResultVector results2;
    if( isTwoReadsOverlap(mergedseq))
		return 1;
    
	//BFS search from 1st to 2nd read via FM-index walk
    while(!m_leaves.empty() && m_leaves.size() <= m_MaxLeaves && m_currentLength <=m_MaxLength)
    {
        // ACGT-extend the leaf nodes via updating existing SA interval
         extendLeaves();
/*              for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
            {
               std::cout <<m_leaves.size()<<" "<< (*iter)->getFullString() <<" Pos="<<m_currentLength<<" :: fwd:"<<(*iter)->fwdInterval.size()<<", rvc:"<<(*iter)->rvcInterval.size()<<", kmersize:"<<m_currentKmerSize<<", kmercount:"<<(*iter)->getKmerCount()<< std::endl;
               // std::cout <<m_leaves.size()<<" "<< (*iter)->getFullString().substr(m_currentLength-m_currentKmerSize) <<" Pos="<<m_currentLength<<" :: fwd:"<<(*iter)->fwdInterval.size()<<", rvc:"<<(*iter)->rvcInterval.size()<<", kmersize:"<<m_currentKmerSize<< std::endl;
            }  */
        // if(m_leaves.size()==1)
        // {
            // mergedseq=(*m_leaves.front()).getFullString();
        // }    
            
            
            //print freq
            // size_t SumOfFreq=0;
            // for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
            // {
                // SumOfFreq = SumOfFreq + (*iter)->fwdInterval.size() + (*iter)->rvcInterval.size(); //accumulate the fwd and rev freq of each extension
            // }
            // std::cout<<m_currentLength<<" "<<SumOfFreq<<" "<<m_currentKmerSize<<" "<<m_leaves.size()<<std::endl;
            // std::cout<<m_currentLength<<" "<<SumOfFreq<<std::endl;
            
            //print pos, kmer, freq  for each branch
            // for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
            // {
                // std::cout<<m_currentLength<<" "<<m_currentKmerSize<<" "<<(*iter)->fwdInterval.size() + (*iter)->rvcInterval.size()<<std::endl;
            // }
            
		if(m_leaves.size()>m_maxUsedLeaves)
            m_maxUsedLeaves=m_leaves.size();
        if(m_leaves.size()>m_MaxLeaves)
        {   
/*             testingOneEND(results2);
            std::string tmp;
            for (size_t i = 0 ; i < results2.size() ;i++)
            {
                //bug fix: m_secondread may be shorter than m_minOverlap
                if(m_secondread.length()>m_minOverlap)
                    tmp=results2[i].thread+m_secondread.substr(m_minOverlap);
                else
                    tmp=results2[i].thread;
                
                size_t cov = calculateKmerCoverage (tmp, m_minOverlap, m_indices.pBWT);
                // size_t cov=results[i].SAICoverage;
                if (  cov > m_maxKmerCoverage )
                {
                    mergedseq=tmp;
                    m_maxKmerCoverage=cov;
                }			
            } */
            
            for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
            {
                size_t cov = (*iter)->getKmerCount();
                if(cov > m_maxKmerCoverage)
                {
                   mergedseq=(*iter)->getFullString();
                   m_maxKmerCoverage=cov;
                }                
            }
        }
		// std::cout << m_currentKmerSize << ":" << m_currentLength << ":" << m_leaves.size() << "\n";	

		//see if terminating string is reached
		if(isTerminated(results))
			break;
    }
		
	//find the path with maximum kmer coverage
	if( results.size()>0 )
	{
		//if multiple paths are bubbles collapsing all together at terminal
		if(results.size()==m_leaves.size()) 
		{
			// std::cout << m_maxUsedLeaves << "\t" << results.size() <<"\n" << results[0].thread << "\n";
			m_isBubbleCollapsed=true;
		}
		std::string tmpseq;
/* 		for (size_t i = 0 ; i < results.size() ;i++)
        {
			//bug fix: m_secondread may be shorter than m_minOverlap
			if(m_secondread.length()>m_minOverlap)
				tmpseq=results[i].thread+m_secondread.substr(m_minOverlap);
			else
				tmpseq=results[i].thread;
			
			size_t cov = calculateKmerCoverage (tmpseq, m_minOverlap, m_indices.pBWT);
			// size_t cov=results[i].SAICoverage;
			if (  cov > m_maxKmerCoverage )
			{
				mergedseq=tmpseq;
				m_maxKmerCoverage=cov;
			}			
		} */
        for (size_t i = 0 ; i < results.size() ;i++)
        {
			size_t cov = results[i].SAICoverage;
			// size_t cov=results[i].SAICoverage;
			if (  cov > m_maxKmerCoverage )
			{
                if(m_secondread.length()>m_minOverlap)
                    tmpseq=results[i].thread+m_secondread.substr(m_minOverlap);
                else
                    tmpseq=results[i].thread;
				mergedseq=tmpseq;
				m_maxKmerCoverage=cov;
			}			
		}
		// for (size_t i=0; i<=mergedseq.length()-m_minOverlap;i+=1)
            // std::cout<<i+1<<" "<<BWTAlgorithms::countSequenceOccurrences(mergedseq.substr(i,m_minOverlap) , m_indices.pBWT )<<std::endl;		
		// std::cout << ">\n" << *m_pQuery << "\n>\n" << reverseComplement(m_secondread);
		// std::cout << mergedseq.length() << "\t" << m_maxKmerCoverage <<  "\t" << (double)m_maxKmerCoverage/mergedseq.length()  << "\n";
		return 1;
    }

	// if(m_leaves.size() > 512){
		// std::cout << m_leaves.size() << "\t" << m_currentLength<< "\t" << m_currentLength-m_pQuery->length()<< "\n";
		// printAll();
		// getchar();
	// }
	
    //Did not reach the terminal kmer
    if(m_leaves.empty())
        return -1;	//high error
    else if(m_currentLength>m_MaxLength)
        return -2;	//exceed search depth
    else if(m_leaves.size() > m_MaxLeaves)
        return -3;	//too much repeats
	else
		return -4;
}
//repeat case -merging round2
int MPSAIntervalTree::mergeRepeat(std::string &mergedseq)
{
    MPSAIntervalNodeResultVector results;
    if( isTwoReadsOverlap(mergedseq))
		return 1;
        
    while(!m_leaves.empty() && m_leaves.size() <= m_MaxLeaves && m_currentLength <=m_MaxLength)
    {
        // ACGT-extend the leaf nodes via updating existing SA interval
        extendLeaves();    
        
        if(m_leaves.size()>m_maxUsedLeaves)
            m_maxUsedLeaves=m_leaves.size();
        
        if(isTerminated(results))
			break;
    }
    
    if( results.size()>0 )
	{
		//if multiple paths are bubbles collapsing all together at terminal
		if(results.size()==m_leaves.size()) 
		{
			// std::cout << m_maxUsedLeaves << "\t" << results.size() <<"\n" << results[0].thread << "\n";
			m_isBubbleCollapsed=true;
		}
		std::string tmpseq;
/* 		for (size_t i = 0 ; i < results.size() ;i++){
			//bug fix: m_secondread may be shorter than m_minOverlap
			if(m_secondread.length()>m_minOverlap)
				tmpseq=results[i].thread+m_secondread.substr(m_minOverlap);
			else
				tmpseq=results[i].thread;
			
			size_t cov = calculateKmerCoverage (tmpseq, m_minOverlap, m_indices.pBWT);
			// size_t cov=results[i].SAICoverage;
			if (  cov > m_maxKmerCoverage )
			{
				mergedseq=tmpseq;
				m_maxKmerCoverage=cov;
			}			
		} */
        for (size_t i = 0 ; i < results.size() ;i++)
        {
			size_t cov = results[i].SAICoverage;
			// size_t cov=results[i].SAICoverage;
			if (  cov > m_maxKmerCoverage )
			{
                if(m_secondread.length()>m_minOverlap)
                    tmpseq=results[i].thread+m_secondread.substr(m_minOverlap);
                else
                    tmpseq=results[i].thread;
				mergedseq=tmpseq;
				m_maxKmerCoverage=cov;
			}			
		}

		return 1;
    }

    //Did not reach the terminal kmer
    if(m_leaves.empty())
        return -1;	//high error
    else if(m_currentLength>m_MaxLength)
        return -2;	//exceed search depth
    else if(m_leaves.size() > m_MaxLeaves)
        return -3;	//too much repeats
	else
		return -4;
}

//On success return the length of merged string
void MPSAIntervalTree::testingOneEND(MPSAIntervalNodeResultVector& results){
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
            std::string STNodeStr = (*iter)->getFullString();
            MPSAIntervalNodeResult STresult;
            STresult.thread=STNodeStr;
			STresult.SAICoverage=(*iter)->getKmerCount();

            //compute the merged pos right next to the kmer on 2nd read.
            results.push_back(STresult);
    }
}

// Print the string represented by every node
void MPSAIntervalTree::printAll()
{
    std::cout << "Print all: \n";
    m_pRootNode->printAllStrings("");
}
void MPSAIntervalTree::filterLowCoverageInterval(STNodePtrList &newLeaves)
{
    size_t maxInterval=0,eachLeaveInterval=0;
    STNodePtrList filterLeaves;// reserve leaves that passed threshold
    for(STNodePtrList::iterator iter = newLeaves.begin(); iter != newLeaves.end(); ++iter)
    {
        eachLeaveInterval = (*iter)->fwdInterval.size()+(*iter)->rvcInterval.size();
        maxInterval = eachLeaveInterval > maxInterval ? eachLeaveInterval : maxInterval;
        if(maxInterval>m_repeatThreshold){
            maxInterval=m_repeatThreshold;
            break;
        }
    }

    for(STNodePtrList::iterator iter = newLeaves.begin(); iter != newLeaves.end(); ++iter)
    {
        eachLeaveInterval = (*iter)->fwdInterval.size()+(*iter)->rvcInterval.size();
        if(maxInterval/10 < eachLeaveInterval)
        {
            filterLeaves.push_back(*iter);
        }
            
    }
    newLeaves.clear();
    newLeaves=filterLeaves;
    
}
// Extend each leaf node
void MPSAIntervalTree::attempToExtend(STNodePtrList &newLeaves)
{
    initCurrentExtendAllFreq();//20160614 chaohung103
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
        extensions = getFMIndexExtensions(*iter);

        // Either extend the current node or branch it
        // If no extension, do nothing and this node
        // is no longer considered a leaf
        if(extensions.size() == 1)
        {
            // Single extension, do not branch
            (*iter)->extend(extensions.front().first);
            (*iter)->fwdInterval=extensions.front().second.interval[0];
            (*iter)->rvcInterval=extensions.front().second.interval[1];
			if((*iter)->fwdInterval.isValid())
            {
                (*iter)->addKmerCount( (*iter)->fwdInterval.size());
                addCurrentExtendAllFreq((*iter)->fwdInterval.size());
            }
			if((*iter)->rvcInterval.isValid())
            { 
                (*iter)->addKmerCount( (*iter)->rvcInterval.size());
                addCurrentExtendAllFreq((*iter)->rvcInterval.size());
            }
			
            newLeaves.push_back(*iter);
        }
        else if(extensions.size() > 1)
        {
            // Branch
            for(size_t i = 0; i < extensions.size(); ++i)
            {
                MPSAIntervalNode* pChildNode = (*iter)->createChild(extensions[i].first);
                pChildNode->fwdInterval=extensions[i].second.interval[0];
                pChildNode->rvcInterval=extensions[i].second.interval[1];
				//inherit accumulated kmerCount from parent
				pChildNode->addKmerCount( (*iter)->getKmerCount() );
				if(pChildNode->fwdInterval.isValid())
                {
                    pChildNode->addKmerCount( pChildNode->fwdInterval.size());
                    addCurrentExtendAllFreq((*iter)->fwdInterval.size());
                }
				if(pChildNode->rvcInterval.isValid()){
                    pChildNode->addKmerCount( pChildNode->rvcInterval.size());
                    addCurrentExtendAllFreq((*iter)->rvcInterval.size());
                }

                newLeaves.push_back(pChildNode);
            }
        }
    }	
}


void MPSAIntervalTree::extendLeaves()
{
    STNodePtrList newLeaves;
	
	//attempt to extend one base for each leave
    attempToExtend(newLeaves);
	//全部葉子加總
    size_t eachLeaveInterval;
    //shrink the SAIntervals in case overlap is larger than read length, which lead to empty newLeaves
    if(!m_kmerMode  &&  newLeaves.empty() )
    {
    
        //refineSAInterval(61);
        eachLeaveInterval=getCurrentExtendAllFreq();
        if(eachLeaveInterval<10)
        {
            for(size_t kmer=m_minOverlap-2;kmer>21 && newLeaves.empty();kmer=kmer-2){
                refineSAInterval(kmer);
                attempToExtend(newLeaves);
            }
        }
        else
        { 
            refineSAInterval(m_minOverlap);
            attempToExtend(newLeaves);
        }
    }	
	
	//extension succeed
    if(!newLeaves.empty()){
		m_currentKmerSize++;
        m_currentLength++;  
	}
    filterLowCoverageInterval(newLeaves);
    m_leaves.clear();
    m_leaves = newLeaves;

	if(!m_leaves.empty() && (m_kmerMode || m_currentKmerSize >= m_maxOverlap) )
    {
        refineSAInterval(m_minOverlap);
        eachLeaveInterval=getCurrentExtendAllFreq();
        if(eachLeaveInterval>m_repeatThreshold * m_leaves.size()){
            refineSAInterval(81);
        }
        // else
            // refineSAInterval(m_minOverlap);
    }
		//refineSAInterval(m_minOverlap);

}

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool MPSAIntervalTree::isTerminated(MPSAIntervalNodeResultVector& results)
{
	bool found = false;

    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        BWTInterval currfwd=(*iter)->fwdInterval;
        BWTInterval currrvc=(*iter)->rvcInterval;

        // assert(currfwd.isValid() || currrvc.isValid());

		//If terminating kmer is a substr, the current SA interval is a sub-interval of the terminating interval
        bool isFwdTerminated = currfwd.isValid() && currfwd.lower >= m_fwdTerminatedInterval.lower
                            && currfwd.upper <= m_fwdTerminatedInterval.upper;
														
        bool isRvcTerminated=currrvc.isValid() && currrvc.lower >= m_rvcTerminatedInterval.lower
                            && currrvc.upper <= m_rvcTerminatedInterval.upper;

        if(isFwdTerminated || isRvcTerminated)
        {
            std::string STNodeStr = (*iter)->getFullString();
            MPSAIntervalNodeResult STresult;
            STresult.thread=STNodeStr;
			STresult.SAICoverage=(*iter)->getKmerCount();

            //compute the merged pos right next to the kmer on 2nd read.
            results.push_back(STresult);
            found =  true;
        }
    }

    return found;
}

bool MPSAIntervalTree::isTwoReadsOverlap(std::string & mergedseq)
{
    //case 1: 1st read sense overlap to 2nd read at exact m_minOverlap bases
    if(BWTInterval::equal(m_pRootNode->fwdInterval, m_fwdTerminatedInterval))
    {
        mergedseq= (*m_pQuery)+m_secondread.substr(m_minOverlap);
        return true;
    }

    //case 2: 1st read sense overlap 2nd read
    std::string secondLeftKmer=m_secondread.substr(0,m_minOverlap);
	//assume overlap can't exceed 100 bp
    size_t pos=m_pQuery->find(secondLeftKmer, m_pQuery->length()>=200?m_pQuery->length()-200:0);
    if(pos!=std::string::npos)	
    {
		//make sure entire suffix of 1st read after pos matches the prefix of 2nd read
		if( m_pQuery->substr(pos) == m_secondread.substr(0, m_pQuery->length()-pos) )
		{
			mergedseq=m_pQuery->substr(0,pos)+m_secondread;
			return true;
		}
    }

    //case 3: 1st read antisense overlap with 2nd read, or 1st read is substr of 2nd read
	//This is rare case and we don't do this in m_kmerMode during island joint
	if(m_kmerMode) return false;
    std::string firstLeftKmer=m_pQuery->substr(0,m_minOverlap);
    pos=m_secondread.find(firstLeftKmer);
	//assume antisense overlap can't exceed 50bp due to rare cases
    if(pos!=std::string::npos && pos <=50)
    {
        //make sure entire suffix of 2nd read after pos matches the prefix of 1st read
		if( m_secondread.substr(pos) ==  m_pQuery->substr(0, m_secondread.length()-pos))
		{
			//return overlapped portion
			mergedseq=m_secondread.substr(pos);
			return true;
		}
    }

    return false;

}

//update SA intervals of each leaf, which corresponds to one-base extension
std::vector<std::pair<std::string, BWTIntervalPair> > MPSAIntervalTree::getFMIndexExtensions(MPSAIntervalNode* pNode)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;

    for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
    {
        char b = BWT_ALPHABET::getChar(i);

        //update forward Interval using extension b
        BWTInterval fwdProbe=pNode->fwdInterval;
        if(fwdProbe.isValid())
            BWTAlgorithms::updateInterval(fwdProbe, b, m_indices.pRBWT);

        //update reverse complement Interval using extension rcb
        BWTInterval rvcProbe=pNode->rvcInterval;
		char rcb=BWT_ALPHABET::getChar(5-i); //T,G,C,A
        if(rvcProbe.isValid())
            BWTAlgorithms::updateInterval(rvcProbe, rcb, m_indices.pBWT);

        size_t bcount = 0;
        if(fwdProbe.isValid())
            bcount += fwdProbe.size();
        if(rvcProbe.isValid())
            bcount += rvcProbe.size();
			
		//min freq at fwd and rvc bwt
        if(bcount >= m_min_SA_threshold)
        {
			// if(bcount>50)
				// std::cout << m_currentKmerSize << ":" << bcount <<"\n";
            // extend to b
            std::string tmp;
            tmp.append(1,b);
            BWTIntervalPair bip;
            bip.interval[0]=fwdProbe;
            bip.interval[1]=rvcProbe;
            out.push_back(std::make_pair(tmp, bip));
        }
    }// end of ACGT

    return out;
}

size_t MPSAIntervalTree::calculateKmerCoverage (const std::string & seq , size_t kmerLength , const BWT* pBWT)
{
	if (seq.length() < kmerLength) return 0;

	size_t cov = 0 ;
	for (size_t i=0; i<=seq.length()-kmerLength;i+=kmerLength/2)
		cov += BWTAlgorithms::countSequenceOccurrences(seq.substr(i,kmerLength) , pBWT );
	
	return cov;
}

// replace each kmer with highest one at each locus
bool MPSAIntervalTree::replaceLowFreqKmer (std::string & seq , size_t kmerLength)
{
	bool changed = false;
	
	for (size_t i=0; i <=seq.length()-kmerLength; i++)
	{
		//Forward kmer should be computed reversely using pRBWT for forward extension
		BWTInterval fwdProbe=BWTAlgorithms::findInterval(m_indices.pRBWT, reverse(seq.substr(i, kmerLength-1)));
		BWTInterval rvcProbe=BWTAlgorithms::findInterval(m_indices.pBWT, reverseComplement(seq.substr(i, kmerLength-1)));
		
		size_t maxcov=0;
		for(int j = 1; j < BWT_ALPHABET::size; ++j) //j=A,C,G,T
		{
			char b = BWT_ALPHABET::getChar(j);

			//update forward Interval using extension b
			if(fwdProbe.isValid())
				BWTAlgorithms::updateInterval(fwdProbe, b, m_indices.pRBWT);

			//update reverse complement Interval using extension rcb
			char rcb=BWT_ALPHABET::getChar(5-i); //T,G,C,A
			if(rvcProbe.isValid())
				BWTAlgorithms::updateInterval(rvcProbe, rcb, m_indices.pBWT);

			size_t bcount = 0;
			if(fwdProbe.isValid())
				bcount += fwdProbe.size();
			if(rvcProbe.isValid())
				bcount += rvcProbe.size();

			if(bcount > maxcov) {
				maxcov = bcount;
				seq.replace(i+kmerLength-1, 1, 1, b);
				changed = true;
			}
		}
	}
	
	return changed;
}

// Refine SA intervals of each leave with a new kmer
void MPSAIntervalTree::refineSAInterval(size_t newKmerSize)
{
	assert(m_currentLength >= newKmerSize);
    initCurrentExtendAllFreq();//20160614 chaohung103
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        // reset the SA intervals using original m_minOverlap
        std::string pkmer = (*iter)->getSuffix(newKmerSize);
		(*iter)->fwdInterval=BWTAlgorithms::findInterval(m_indices.pRBWT, reverse(pkmer));
		(*iter)->rvcInterval=BWTAlgorithms::findInterval(m_indices.pBWT, reverseComplement(pkmer));
        
        addCurrentExtendAllFreq((*iter)->fwdInterval.size());
        addCurrentExtendAllFreq((*iter)->rvcInterval.size());
    }

	m_currentKmerSize=newKmerSize;
}

/***Dead code***/

// Remove leaves with two or more same kmers
void MPSAIntervalTree::removeLeavesByRepeatKmer()
{
    STNodePtrList newLeaves;    
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        /*
        GAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTG
        TGGATTCCAGATTGTTCGAGGAGAATTTGGTGGAGCTACGCGGGATCGAACCGCGGACCTCTTGCATGCCATGCAAGCGCTCTCCCAGCTGAGCTATAACCCCTTGGATTCCAGATTGTTCGAGGAGAATTTGGTGGAGCTACGCGGGATCGAACCGCGGACCTCTTGCATGCCATGCAAGCGCTCTCCCAGCTGAGCTATAACC
        GAGAGGGACTCGAACCCTCACACCCGGGGGGCACTAACACCTGAAGCTAGCGCGTCTACCAATTCCGCCACCTTCGCACATCGGGTTATCAGTCTGGATTTACATGCTGTCTGATAAAAGCATGGTGCGAAGAGAGGGACTCGAACCCTCACACCCGGGGGGCACTAACACCTGAAGCTAGCGCGTCTACCAATTCCGCCACCTTCGCACATCGGGTTATCAGTCTGGATTT
        GCATATCCATCCCACCAGCACATCGACCTATCGACTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATC
        CATCGGCGTCAGCCTGCTGGGCTTCACCCATCAGGGCAACAAGTGGCTGTGGCAGCAGGCCAGGGCCGCTCTTCCCTCCCTCAAGGGGGAGCTGGTGGCGGGGGGGGGGGGGGGGGGGGGGGGGGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
        */
        std::string STNodeStr = (*iter)->getFullString();
        std::string fwdrepeatunit = STNodeStr.substr(STNodeStr.size()-m_minOverlap);
        std::string revrepeatunit = reverseComplement(fwdrepeatunit);
        size_t index1=STNodeStr.find(fwdrepeatunit);
        size_t index2=STNodeStr.find(revrepeatunit);

        if(index1 == (STNodeStr.size()- m_minOverlap) && index2 == std::string::npos)
        {
            newLeaves.push_back(*iter);
        }
    }

    m_leaves=newLeaves;
}
