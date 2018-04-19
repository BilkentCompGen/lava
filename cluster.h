#ifndef UNTITLED_CLUSTER_H
#define UNTITLED_CLUSTER_H

#include "stringOp.h"
#include "read.h"
#include "sv.h"

using namespace std;

class Cluster {
public:
    int id;
    int32_t chromosome;
    vector<SV *> svs, crit, rCrit;
    int cstart1, cstart2, cend1, cend2, tol;
    int num_of_edges;

    Cluster () {
	id = -1;
	chromosome = -1;
	cstart1 = 0;
	cstart2 = 0;
	cend1 = 0;
	cend2 = 0;
	tol = 0;
	num_of_edges = 0;
    }
};

extern double quasi_clique_lambda, quasi_clique_gamma;
extern int minClusterSize, minSupport;
extern double depthThreshold, ovp;
extern int filterSegDups, filterSats, checkReps;

void addCluster(vector<Cluster> *clusters);

/*
 * Updates break points of the cluster sv.
 */
int updateBpAsCentroid( Cluster *cl, SV *sv, int updateType );
int updateBpAsOutermost( Cluster *cl, SV *sv, int updateType );

/*
 * Updates the number of edges of a cluster and degrees of its SVs according to the
 * newly added or removed sv. Returns 1 if update is successful, 0 otherwise.
 */
int updateCluster( Cluster *cl, SV *sv);

/*
 * First, checks whether the chromosome of the new SV is the same as the cluster's.
 * Then, checks whether the new SV has an edge to all the nodes in the critical
 * addition set of the cluster. If it does, adds the new SV and returns 1. If not,
 * returns 0. Returns -1 in case of errors.
 */
int addSV( Cluster *cl, SV *newSV );

/*
 * Removes an SV from a cluster.
 */
int removeSV( Cluster *cl, SV *sv);

/*
 * Clusters the PacBio SVs using a quasi clique approximation algorithm.
 */
int clusterSVs(vector<SV> *svs, vector<Cluster> *clusters);

void filterDups( vector<SV> *svsTemp, vector<int> sizeArr, vector<int> supportArr );
int filterSVs( vector<SV> *svs, vector<Cluster> *clusters, vector<SV> *finalSVs, const vector<short> *dV, const double meanReadDepth, const int currentChr, ofstream *out );

#endif //UNTITLED_CLUSTER_H
