#ifndef UNTITLED_SV_H
#define UNTITLED_SV_H

#include "read.h"
#include "stringOp.h"

using namespace std;

/*** Minimum and maximum size of deletions to be detected ***/
#define MIN_DEL_LEN 10000
#define MAX_DEL_LEN 10000000

/*** Minimum and maximum size of the inversions to be detected ***/
#define MIN_INV_LEN 50000
#define MAX_INV_LEN 10000000

/*** Quasi-clique lambda and gamma parameters for deletions ***/
#define DEL_LAMBDA 0.5
#define DEL_GAMMA 0.6

/*** Quasi-clique lambda and gamma parameters for inversions ***/
#define INV_LAMBDA 0.2
#define INV_GAMMA 0.3

/*** Desired percentage of overlap between two deletions - Determines whether there is an edge between two vertices ***/
#define DEL_CLUSTER_OVERLAP_PERCENTAGE 0.50
#define DEL_VERIFICATION_OVERLAP_PERCENTAGE 0.70

/*** Desired percentage of overlap between two inversions - Determines whether there is an edge between two vertices ***/
#define INV_CLUSTER_OVERLAP_PERCENTAGE 0.50
#define INV_VERIFICATION_OVERLAP_PERCENTAGE 0.95

/*** Clusters with smaller number of deletions are filtered out ***/
#define DEL_MIN_CLUSTER_SIZE 3

/*** Clusters with smaller number of inversions are filtered out ***/
#define INV_MIN_CLUSTER_SIZE 20

/*** Clusters with smaller number of Illumina matches are filtered out ***/
#define DEL_MIN_SUPPORT 15

/*** Clusters with smaller number of Illumina matches are filtered out ***/
#define INV_MIN_SUPPORT 70

/*** Deletions with higher read depth are filtered out ***/
#define DEL_DEPTH_THRESHOLD 2

/*** Inversions with higher read depth are filtered out ***/
#define INV_DEPTH_THRESHOLD 12

/*** Min quality score to consider a read ***/
#define PACBIO_MIN_QUAL 0
#define ILMN_MIN_QUAL 0

/*** Decide whether to check these regions for deletions or not ***/
#define DEL_FILTER_SEGMENTAL_DUPLICATIONS 1
#define DEL_FILTER_SATELLITES 1

/*** Decide whether to check these regions for inversions or not ***/
#define INV_FILTER_SEGMENTAL_DUPLICATIONS 0
#define INV_FILTER_SATELLITES 1

/*** Prints repeat type and location within each detected deletion if this is 1 ***/
#define DEL_CHECK_REPEATS 0

/*** Prints repeat type and location within each detected inversion if this is 1 ***/
#define INV_CHECK_REPEATS 0

/*** Parameters for detecting deletions within the cigar field of reads ***/
#define INTERFERING_BPS 10
#define MIN_DEL_STARTING_LEN 1

typedef unordered_multimap<std::string, Read> Unmmap;
typedef allocator<pair<const std::string, Read> > alloc;

typedef unordered_multimap<std::string, vector<vector<int> > > Altmap;
typedef allocator<pair<const std::string, vector<vector<int> > > > altAlloc;

class Params {
private:
    Params()
    {
	pacbio_bam_file = "";
	ilmn_bam_file = "";
	name_of_reference = "";
	output = "";
	min_del_len = MIN_DEL_LEN;
	min_inv_len = MIN_INV_LEN;
	min_del_cluster_size = DEL_MIN_CLUSTER_SIZE;
	min_del_support = DEL_MIN_SUPPORT;
	min_inv_cluster_size = INV_MIN_CLUSTER_SIZE;
	min_inv_support = INV_MIN_SUPPORT;
	pacbio_min_qual = PACBIO_MIN_QUAL;
	ilmn_min_qual = ILMN_MIN_QUAL;
    }

    static Params *instance;

public:
    static Params* getInstance()
    {
	if (instance == NULL) {
	    instance = new Params;
	}

	return instance;
    }

    /*** required arguments ***/
    string pacbio_bam_file;
    string ilmn_bam_file;
    string name_of_reference;

    /*** optional arguments ***/
    string output;

    int min_del_len;
    int min_inv_len;

    int min_del_cluster_size;
    int min_del_support;
    int min_inv_cluster_size;
    int min_inv_support;

    int pacbio_min_qual;
    int ilmn_min_qual;
};

struct position {
    uint64_t start1;
    uint64_t start2;
    uint64_t end1;
    uint64_t end2;
};

enum svtype {
    none,
    deletion,
    inversion
};

class SV {
  public:
    svtype type;
    int32_t chromosome;
    struct position pos;
    int tolerance;
    int clusterID;
    int clusterIndex;
    int degree;
    int vcfID;
    string readID;

    SV () {
	type = none;
	chromosome = -1;
	pos.start1 = 0;
	pos.start2 = 0;
	pos.end1 = 0;
	pos.end2 = 0;
	tolerance = 0;
	clusterID = -1;
	clusterIndex = -1;
	degree = 0;
	vcfID = -1;
	readID = " ";
    }
};

int validateRead( Read *rd, int currentChr, sonic *snc, int type );

int percOvlp( const int start1, const int start2, const int end1, const int end2, double overlapPcg );
int toleranceOvlp( const SV *s1, const SV *s2 );

int delOvlp( const SV *s1, const SV *s2, double overlapPcg );
int invOvlp( const SV *s1, const SV *s2, double overlapPcg );

int assignEdge( const SV *s1, const SV *s2, double overlapPcg);

int mergeInv( const SV *s1, const SV *s2);

int delCigar( const Read *rd, uint64_t *start, uint64_t *len, uint64_t *tol );

void detectLowDepth( vector<SV> *deletions, const vector<SV> *il_deletion, const int currentChr, sonic *snc, const vector<short> *depthVector, ofstream *out );

int verifyInversions( const vector<SV> *forwardSplits, const vector<SV> *reverseSplits, vector<SV> *inversions, double ovpc, int merge );

int findIlluminaDel( vector<SV> *deletions, const bam1_t *aln, const int currentChr, sonic *snc );
int findIlluminaInv( vector<SV> *inversions, vector<SV> *forwardPairs, vector<SV> *reversePairs, const bam1_t *aln, const int currentChr, sonic *snc );

int findPacbioDel ( vector<SV> *deletions, const Read *cur_rd, const Read *next_rd, const int currentChr, sonic *snc );
int findPacbioInv ( vector<SV> *forwardSplits, vector<SV> *reverseSplits, const Read *cur_rd, const Read *next_rd, const int currentChr, sonic *snc );
int findPacbioSVs( Unmmap *candidates, vector<SV> *deletions, vector<SV> *inversions, const int currentChr, sonic *snc );

#endif //UNTITLED_SV_H
