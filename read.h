#ifndef UNTITLED_READ_H
#define UNTITLED_READ_H

#include "stringOp.h"

using namespace std;

class Read {
public:
    int32_t rname;
    int32_t pos;
    uint16_t flag;
    int32_t m_rname;
    int32_t m_pos;
    int32_t ncigar;
    uint32_t *cigar;
    string qname;

    Read () {
	rname = -1;
	pos = 0;
	flag = 0;
	m_rname = -1;
	m_pos = 0;
	ncigar = 0;
	qname = " ";
    }
};

vector<int> totalOpCount( const Read *rd, int numofargs, ... );

/*
 * Finds the length (in terms of number of bases) of the soft clipped (sc) regions at the beginning and end of a read.
 * Outputs a pair of integers, first one denoting the length of the sc at the beginning of the read and second one
 * denoting the length of the sc at the end of the read.
 */
void softClipLengths( const Read *rd, vector<int> *sclengths );

/*
 * Computes read reference span (#M + #D)
 */
int referenceSpan( const Read *rd );

/*
 * Computes read length (#M + #I)
 */
int readLength( const Read *rd );

/*
 * Computes illumina mate reference span
 */
int mateReferenceSpan( const bam1_t *aln );

/*
 * Computes alternative reference span (#M + #D)
 */
int altRefSpan(const char *altCigar);

/*
 * Initializes variables for read depth calculation
 */
int _readDepth_init( long long unsigned int *totalReadLengths, vector<short> *depthVector, const int currentChr, const bam_hdr_t *bamHdr );

/*
 * Adds the number of all matching bases of the read with the reference to the depth vector.
 * Also increments the global read lengths by number of matches. Returns number of matches
 * on success. If the read is in satellite region, it is not included in depth and
 * addRdToDepth returns -1.
 */
int addRdToDepth( const Read *rd, vector<short> *depthVector, unsigned long long int *totalReadLengths );

/*
 * Finds all alternative mappings of a read. Pushes alternative mappings in the current
 * chromosome to a vector, writes alternative mappings in other chromosomes to a file.
 * Returns 1 on success, 0 otherwise.
 */
int findAltMappings (const int currentChr, const string readID, string xa, vector <vector<string> > *altMappings, ofstream *alts );

#endif //UNTITLED_READ_H
