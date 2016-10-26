///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// MPExtensionProcess.cpp - Implementation of FM-index walk and kmerization of PE reads
//
#include "MPExtensionProcess.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include <iomanip>
#include "MPSAIntervalTree.h"
using namespace std;

//#define KMER_TESTING 1


MPExtensionProcess::MPExtensionProcess(const MPExtensionParameters params) : m_params(params)
{
}

MPExtensionProcess::~MPExtensionProcess()
{

}




//trimmatepair by chaohung 20151017
MPExtensionResult MPExtensionProcess::MatepairWorkflow(const SequenceWorkItemPair& workItemPair)
{
	MPExtensionResult result;
	std::string mergedseq1,mergedseq2;
	int longRead=0;
    const size_t RepeatKmerFreq = m_params.kd.getMedian()*1.3;
    //const size_t FreqThreshold = m_params.kd.findFirstLocalMinimum();
	//get parameters123
	size_t kmerLength = m_params.kmerLength ;
	size_t threshold = (size_t)CorrectionThresholds::Instance().getRequiredSupport(0)-1;
    
	std::string seqFirst  = workItemPair.first.read.seq.toString() ;
	std::string seqSecond = workItemPair.second.read.seq.toString();
	
	//Trim head and tail of both ends containing errors
    seqFirst = trimRead(seqFirst, kmerLength,threshold,m_params.indices,&result.head1, m_params.FreqThreshold);
	seqSecond = trimRead(seqSecond, kmerLength,threshold,m_params.indices,&result.head2, m_params.FreqThreshold);
	
	if(!seqFirst.empty() && !seqSecond.empty())//both 1st and 2nd not empty
	{
		result.tryMerge=true;
		std::string firstKRstr = reverseComplement(seqFirst.substr(0));
		std::string secondKRstr  = seqSecond.substr(0);

		
		// maxOverlap is limited to 90% of read length which aims to prevent over-greedy search
        size_t maxOverlap = m_params.maxOverlap!=-1?m_params.maxOverlap:
											((workItemPair.first.read.seq.length()+workItemPair.second.read.seq.length())/2)*0.95;

	
		// Walk from the 1st end to 2nd end
		MPSAIntervalTree SAITree1(&firstKRstr, m_params.minOverlap, maxOverlap, m_params.maxInsertSize, m_params.maxLeaves,
											m_params.indices, RepeatKmerFreq, secondKRstr);
			
		longRead=SAITree1.mergeTwoReads(mergedseq1);
        
             
  		 if(longRead==-3)//round 2 - merging repeat case
         {  //第一次走的路線要rvc 才能當作第二次的終點
            
            firstKRstr = reverseComplement(mergedseq1.substr(0,mergedseq1.length()-2));
            //firstKRstr = reverseComplement(mergedseq1.substr(0,mergedseq1.length()-m_params.minOverlap));
            secondKRstr = reverseComplement(secondKRstr);
            MPSAIntervalTree SAITree2(&secondKRstr, m_params.minOverlap, maxOverlap, m_params.maxInsertSize-mergedseq1.length()+30, m_params.maxLeaves,
                                                m_params.indices, RepeatKmerFreq, firstKRstr);
                  
            longRead=SAITree2.mergeRepeat(mergedseq1);
			
            // std::cout<<"second round"<<std::endl;
        } 
			
		// std::cout<<workItemPair.first.read.id.substr (0, workItemPair.first.read.id.find('/') )<<" : "<<"longRead="<<longRead<<", shortIS="<<shortIS<<std::endl;

        
		if(longRead>-2 && !mergedseq1.empty())
		{
			if(mergedseq1.length()> (unsigned)m_params.minInsertSize){
				result.mergeMP =true;
				result.correctSequence = mergedseq1;
				result.longmergeCount=true;
				return result;
			}
			else
			{
				result.shortmergeCount = true;
                result.whyLongFail=longRead=-4;
				result.correctSequence = mergedseq1;
				return result;
			}
			
		}
		else
		{
			result.trim = true ;
			result.whyLongFail=longRead;
			result.correctSequence = seqFirst ;
			result.correctSequence2 = seqSecond;
			return result;
		}
		
	}
	else
	{
		result.polluted = true;
		result.pollutedSequence = workItemPair.first.read.seq.toString() ;
		result.pollutedSequence2 = workItemPair.second.read.seq.toString();
		return result;
		
	}
    
	return result;
}
//
MPExtensionResult MPExtensionProcess::process(const SequenceWorkItem& workItem)
{
	MPExtensionResult result = correct(workItem);
	return result;
}

MPExtensionResult MPExtensionProcess::correct(const SequenceWorkItem& /*workItem*/)
{

	MPExtensionResult result;
	return result;
}




// return complexity of seq, default: 0.9
bool  MPExtensionProcess::isLowComplexity (std::string seq , float & GCratio)
{
	size_t seqLen = seq.length();
	size_t countG =0 ;
	size_t countC =0 ;
	size_t countT =0 ;
	size_t countA =0 ;

	for (size_t i=0; i<seqLen; i++)
	{
		switch(seq[i]){
			case 'A': countA ++ ;break;
			case 'T': countT ++ ;break;
			case 'C': countC ++ ;break;
			case 'G': countG ++ ;break;
			default:  assert(false);
		}
	}

	GCratio = (float)(countG+countC)/seqLen ;

	if (  ((float) countA/seqLen >=0.9 ) || ((float) countT/seqLen >=0.9 )
	   || ((float) countC/seqLen >=0.9 ) || ((float) countG/seqLen >=0.9 ) )
	   return true;

	return false;

}



//search for a strong interval having high-frequent kmers at both strands 
bool MPExtensionProcess::isIntervalExistStrongKmer (std::pair<size_t,size_t> interval,std::vector<size_t> & countQualified)
{
	for(size_t i = interval.first ; i<=interval.second ; i++ )
	{
		//find a high-frequent kmer 
		if (countQualified.at(i)== 2 ) return true;
	}
	return false;
}

//Determine reliability between two intervals defined by existence of strong kmer at one strand
bool MPExtensionProcess::isPathReliable(std::pair<size_t,size_t> intervalX, std::pair<size_t,size_t> intervalY,std::vector<size_t> & countQualified)
{
	//Two adjacent strong intervals are assumed to be reliable
	if (intervalX.second+1==intervalY.first) return true;

	size_t start = intervalX.second + 1 ;
	size_t end = intervalY.first-1;
	assert (start<=end);

	//if the path exists strong kmer at one strand, it is reliable.
	//Otherwise, 
	for (size_t i =start  ; i <=end ;i++)
		if (countQualified[i]==0) return false;

	return true;
}

//Merge back two splitted intervals if there exists strong kmer link at one strand
bool MPExtensionProcess::isIntervalMerge (std::vector< std::pair<size_t,size_t> > & intervals , std::vector<size_t> & countQualified )
{
	std::vector<bool> stongInterval(intervals.size());
	size_t count = 0 ;
	for (size_t i =0 ;i <intervals.size();i++)
	{
		stongInterval.at(i)= isIntervalExistStrongKmer (intervals[i],countQualified);
		if (stongInterval[i]) count ++ ;
	}
	if (count <2 ) return false;
	else	//This read exists two strong kmers in two intervals, may be over splitted
	{
		int s = -1 ;
		int e = -1 ;
		for (int i=0 ;i<(int)intervals.size();i++)
		{
			//Find the first strong interval
			if (stongInterval.at(i)&& s<0)
			{
				s=i;
				continue;
			}

			//Find the second strong interval
			if (stongInterval[i])
			{
				e=i;
				//check if there exist strong kmers at one strand between them
				if (isPathReliable(intervals[s],intervals[e],countQualified))
				{
					//Extend 1st interval end to 2nd end
					intervals[s].second = intervals[e].second;
					//Erase the intervals after 1st interval (s+1) till 2nd interval e
					intervals.erase(intervals.begin()+s+1,intervals.begin()+e+1);
					return true;
				}
				s=e;
			}
		}
	}
	return false;

}





//trimmatepair by chaohung 20151017
std::string MPExtensionProcess::trimRead( std::string readSeq ,size_t kmerLength , size_t /*threshold*/ ,BWTIndexSet & index,bool *H, const size_t freqThrshold)
{


	
	//trim
	int kmerPos=0,head=0, tail=0, lastKmer=readSeq.length()-kmerLength;
	
	

	

	kmerPos=0;
	for(kmerPos = kmerPos ;kmerPos<=lastKmer;kmerPos++){
		if (BWTAlgorithms::countSequenceOccurrences( readSeq.substr(kmerPos,kmerLength), index)<=freqThrshold)
		{
			continue;
		}
		else
		{
				head=kmerPos;
                *H=true;
				break;
		}
	}
	if(kmerPos<=lastKmer){
		for (kmerPos = kmerPos ; kmerPos <= lastKmer ; kmerPos++ )
		{
			if (BWTAlgorithms::countSequenceOccurrences( readSeq.substr(kmerPos,kmerLength), index)>freqThrshold)
			{
				if(kmerPos==lastKmer){
					tail=kmerPos;
				}	
				else
					continue;
			}
			else
			{								
				tail=kmerPos-1;
				break ;
				
			}
		}
	}
   
    
    
	
	//all kemrs are dirty , return empty
	if (!*H)
		return "";
	else
		return readSeq.substr(head,tail-head+kmerLength);//(head,tail+kmerLength);

}


size_t MPExtensionProcess::numNextKmer(std::string kmer , NextKmerDir dir ,BWTIndexSet & index, size_t threshold)
{
	size_t num = 0 ;
	char nBases[4] = {'A','T','C','G'} ;
	int kmerLength = kmer.length() ;

	for (size_t i = 0 ; i < 4 ; i++)
	{
		std::string next_mer;
		if (dir == NK_START) next_mer = nBases[i] + kmer.substr(0,kmerLength-1);
		else if  (dir == NK_END) next_mer = kmer.substr(1,kmerLength-1) + nBases[i];

		if ( BWTAlgorithms::countSequenceOccurrences(next_mer,index) >= threshold) num++;
	}
	return num ;
}


bool MPExtensionProcess::isSimple (std::string Lkmer, std::string Rkmer, BWTIndexSet & index, size_t threshold)
{
	size_t LKmerPathCount = numNextKmer(Lkmer, NK_END, index, threshold);
	size_t RKmerPathCount = numNextKmer(Rkmer, NK_START, index, threshold);
	if ( LKmerPathCount == 1 &&   RKmerPathCount == 1 )  
		return true;
	else 
		return false;
	
}



//
//
//
MPExtensionPostProcess::MPExtensionPostProcess(std::ostream* pCorrectedWriter,
																				std::ostream* pDiscardWriter,
																				const MPExtensionParameters params) :
																													m_pCorrectedWriter(pCorrectedWriter),
																													m_pDiscardWriter(pDiscardWriter),
																													m_params(params),
																													m_kmerizePassed(0),
																													m_mergePassed(0),
                                                                                                                    m_qcFail(0)                                                                                                                    
																				{
																					//m_ptmpWriter = createWriter("NoPESupport.fa");
																				}

//trimmatepair by chaohung 20151030
MPExtensionPostProcess::MPExtensionPostProcess(std::ostream* pCorrectedWriter,
																				std::ostream* pDiscardWriter,
																				std::ostream* pLongReadWriter,
																				std::ostream* pShortReadWriter,
																				const MPExtensionParameters params) :
																													m_pCorrectedWriter(pCorrectedWriter),
																													m_pDiscardWriter(pDiscardWriter),
																													m_pLongReadWriter(pLongReadWriter),
																													m_pShortReadWriter(pShortReadWriter),
																													m_params(params),
																													m_kmerizePassed(0),
																													m_mergePassed(0),
                                                                                                                    //trimmatepair by chaohung 20151017
                                                                                                                    m_qcFail(0),
                                                                                                                    m_trimPassed(0),
                                                                                                                    m_r1_pass(0),
                                                                                                                    m_r2_pass(0),
																													m_totalMerge(0),
																													m_longSucc(0),
																													m_HighError(0),
																													m_exceedSearchDepth(0),
																													m_Repeat(0),
																													m_Case4(0)
																				{
																					//m_ptmpWriter = createWriter("NoPESupport.fa");
																				}
//
MPExtensionPostProcess::~MPExtensionPostProcess()
{
    //trimmatepair by chaohung 20151017
		/*std::cout << "Reads trimmed from head success: " << s_numHeadsuccess << "\n";//By_Chaohung_2015_01_27
		std::cout << "Reads trimmed from tail success: " << s_numTailsuccess << "\n";//By_Chaohung_2015_01_27
		std::cout << "Reads trimmed from head fail: " << s_numHeadfail << "\n";//By_Chaohung_2015_01_27
		std::cout << "Reads trimmed from tail fail: " << s_numTailfail << "\n";//By_Chaohung_2015_01_27
        */
        


		std::cout << "    Reads total merge: " << m_totalMerge << "\n";
		std::cout << "        -Long Merge      : " << m_longSucc << "\n";
		std::cout << "            -Long merge fail : " <<"\n";
		std::cout << "                    -HighError          : " << m_HighError << "\n";
		std::cout << "                    -Exceed Search Depth: " << m_exceedSearchDepth << "\n";
		std::cout << "                    -Repeat             : " << m_Repeat << "\n";
		std::cout << "                    -Too short(<0.125*MAX_insertSize)    : " << m_Case4 << "\n";
		std::cout << "    Reads failed to trim: " << m_qcFail << "\n";

}


//
void MPExtensionPostProcess::process(const SequenceWorkItem& item, const MPExtensionResult& result)
{
	// Determine if the read should be discarded
	bool readQCPass = true;
	if (result.kmerize)
	{
		m_kmerizePassed += 1;
	}
	else
	{
		readQCPass = false;
		m_qcFail += 1;
	}

	SeqRecord record = item.read;
	record.seq = result.correctSequence;

	if (result.correctSequence.empty() && result.kmerizedReads.empty()) ;

	else if (result.kmerize)
	{
		if (!result.correctSequence.empty())
			record.write(*m_pCorrectedWriter);


		for (size_t i=0 ; i< result.kmerizedReads.size() ; i++)
		{
			record.seq = result.kmerizedReads[i];
			record.writeFasta(*m_pDiscardWriter,i);
		}

	}
	else if  (readQCPass || m_pDiscardWriter == NULL)
	{
		record.write(*m_pCorrectedWriter);
	}
	else
	{
		record.write(*m_pDiscardWriter);
	}

}

// Writting results of FMW_HYBRID and FMW_MERGE
//trimmatepair by chaohung 20151017
void MPExtensionPostProcess::process(const SequenceWorkItemPair& itemPair, const MPExtensionResult& result)
{

    //trimmatepair by chaohung 20151017

        if(result.head1)
            m_r1_pass+=1;
        if(result.head2)
            m_r2_pass+=1;
		if(result.trim)
			m_trimPassed += 2;
		if(result.tryMerge)
			m_totalMerge+=1;
		if(result.longmergeCount)
			m_longSucc+=1;
        else
            switch(result.whyLongFail){
				case -1:m_HighError+=1;break;
				case -2:m_exceedSearchDepth+=1;break;
				case -3:m_Repeat+=1;break;
                case -4:m_Case4+=1;break;
				default:
				break;
			}
		if(result.polluted)
			m_qcFail+=2;


	SeqRecord firstRecord  = itemPair.first.read;
	SeqRecord secondRecord  = itemPair.second.read;

    //trimmatepair by chaohung 20151017

		if(result.mergeMP)
		{
			SeqRecord mergeRecord ;
			if(result.longmergeCount)
			{
				mergeRecord.id = firstRecord.id.substr (0, firstRecord.id.find('/') ) ;
				mergeRecord.seq = result.correctSequence;
				mergeRecord.write(*m_pLongReadWriter);
			}
			if(result.shortmergeCount)
			{
				mergeRecord.id = firstRecord.id.substr (0, firstRecord.id.find('/') ) ;
				mergeRecord.seq = result.correctSequence;
				mergeRecord.write(*m_pShortReadWriter);
			}
		}
		if (result.trim) //result.trim: only merge fail will redirect here. ;result.tryMerge:all trim read will redirect here
		{	
            std::stringstream ss;
			SeqRecord raw;
			SeqRecord raw2;
			SeqRecord trimRecord;
			SeqRecord trimRecord2;

            //matepair pair-end merge failure output raw and afterTrim read to File
                 //raw r1
                ss.str("");
                ss<<result.whyLongFail<<" ";//Long Merge Failure reason
                raw.id = firstRecord.id.substr(0)+ " "+ss.str();
                ss.str("");
                raw.seq= firstRecord.seq.toString();
                raw.write(*m_pCorrectedWriter);
                //trimmed 1       
                trimRecord.id = firstRecord.id.substr(0, firstRecord.id.find('/') )+"/1";// /trimmed Kmer= "+ss.str();
                trimRecord.seq = result.correctSequence;
                trimRecord.write(*m_pCorrectedWriter);
                //origin 2
                ss.str("");
                ss<<result.whyLongFail<<" ";
                raw2.id = secondRecord.id.substr(0)+ " "+ss.str();
                ss.str("");                
                raw2.seq= secondRecord.seq.toString();
                raw2.write(*m_pCorrectedWriter);
                //trimmed 2
                trimRecord2.id = secondRecord.id.substr(0, firstRecord.id.find('/') )+"/2";// /trimmed Kmer= "+ss.str();
                trimRecord2.seq = result.correctSequence2;
                trimRecord2.write(*m_pCorrectedWriter); 
            
		}
		if (result.polluted)
		{
			std::stringstream ss;
			SeqRecord pollutedRecord;
			SeqRecord pollutedRecord2;
			//ss.str("");
			/* for(int i=0;i<result.c1;++i)
				ss<<result.cs1[i]<<" "; */
            pollutedRecord.id = firstRecord.id.substr (0, firstRecord.id.find('/') )+ "/1";// origin frqs: "+ss.str();
			//pollutedRecord.id = firstRecord.id.substr (0, firstRecord.id.find('/') )+ "/1";// /origin frqs:"+ss.str();
			pollutedRecord.seq = result.pollutedSequence;
			pollutedRecord.write(*m_pDiscardWriter);
			//ss.str("");
			/* for(int i=0;i<result.p1;++i)
				ss<<result.ps1[i]<<" "; */
            pollutedRecord2.id = secondRecord.id.substr(0, secondRecord.id.find('/') )+ "/2";// origin frqs: "+ss.str();
			//pollutedRecord2.id = secondRecord.id.substr(0, secondRecord.id.find('/') )+ "/2";// /origin frqs:"+ss.str();
			pollutedRecord2.seq = result.pollutedSequence2;
			pollutedRecord2.write(*m_pDiscardWriter);
		}
	
}
