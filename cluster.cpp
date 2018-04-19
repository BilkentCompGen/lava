
//
// Created by Ezgi Ebren on 21.03.2018.
//

#include "cluster.h"

double quasi_clique_lambda, quasi_clique_gamma;
double depthThreshold, ovp;
int minClusterSize, minSupport;
int filterSegDups, filterSats, checkReps;

void addCluster(vector<Cluster> *clusters)
{
    Cluster cl;

    cl.id = (*clusters).size();

    (*clusters).push_back(cl);
}

int updateBpAsCentroid( Cluster *cl, SV *sv, int updateType )
{
    uint64_t tol;
    int numOfSVs = (*cl).svs.size();

    double additionMultiplier = (double)(numOfSVs-1) / (double)(numOfSVs);
    double removalMultiplier = (double)(numOfSVs+1) / (double)(numOfSVs);

    if (updateType == 0)
    {
	//Update the centeroid of the cluster after addition
        (*cl).cstart1 = round( ( (double)((*cl).cstart1) * additionMultiplier ) + (double)((*sv).pos.start1) / numOfSVs );
        (*cl).cstart2 = round( ( (double)((*cl).cstart2) * additionMultiplier ) + (double)((*sv).pos.start2) / numOfSVs );
        (*cl).cend1 = round( ( (double)((*cl).cend1) * additionMultiplier ) + (double)((*sv).pos.end1) / numOfSVs );
        (*cl).cend2 = round( ( (double)((*cl).cend2) * additionMultiplier ) + (double)((*sv).pos.end2) / numOfSVs );
        (*cl).tol = round( ( (double)((*cl).tol) * additionMultiplier ) + (double)((*sv).tolerance) / numOfSVs );
    }
    else if (updateType == 1)
    {
        //Update the centeroid of the cluster after removal
	(*cl).cstart1 = round( ( (double)((*cl).cstart1) * removalMultiplier ) - (double)((*sv).pos.start1) / numOfSVs );
	(*cl).cstart2 = round( ( (double)((*cl).cstart2) * removalMultiplier ) - (double)((*sv).pos.start2) / numOfSVs );
	(*cl).cend1 = round( ( (double)((*cl).cend1) * removalMultiplier ) - (double)((*sv).pos.end1) / numOfSVs );
	(*cl).cend2 = round( ( (double)((*cl).cend2) * removalMultiplier ) - (double)((*sv).pos.end2) / numOfSVs );
	(*cl).tol = round( ( (double)((*cl).tol) * removalMultiplier ) - (double)((*sv).tolerance) / numOfSVs );
    }
    else {
	return -1;
    }

    return 1;
}

/** When using this, update start locations as INT_MAX **/
int updateBpAsOutermost( Cluster *cl, SV *sv, int updateType )
{
    SV *tempsv;

    int s1, s2, e1, e2;
    uint64_t tol;
    int numOfSVs = (*cl).svs.size();

    if (updateType == 0)
    {
	//Update the start-end of the cluster after addition
        (*cl).cstart1 = ((*sv).pos.start1 < (*cl).cstart1) ? (*sv).pos.start1 : (*cl).cstart1;
        (*cl).cstart2 = ((*sv).pos.start2 < (*cl).cstart2) ? (*sv).pos.start2 : (*cl).cstart2;
        (*cl).cend1 = ((*sv).pos.end1 > (*cl).cend1) ? (*sv).pos.end1 : (*cl).cend1;
        (*cl).cend2 = ((*sv).pos.end2 > (*cl).cend2) ? (*sv).pos.end2 : (*cl).cend2;

        tol = ((*cl).tol * (numOfSVs - 1)) + (*sv).tolerance;
        (*cl).tol = tol / numOfSVs;
    }
    else if (updateType == 1)
    {
	//Update the start-end of the cluster after removal
        s1 = (*sv).pos.start1 == (*cl).cstart1;
        s2 = (*sv).pos.start2 == (*cl).cstart2;
        e1 = (*sv).pos.end1 == (*cl).cend1;
        e2 = (*sv).pos.end2 == (*cl).cend2;
        if (s1 || s2 || e1 || e2)
        {
            for (int ind = 0; ind < (*cl).svs.size(); ind++)
            {
                tempsv = (*cl).svs[ind];
                if (s1) {
                    (*cl).cstart1 = ( (*tempsv).pos.start1 < (*cl).cstart1 ) ? (*tempsv).pos.start1 : (*cl).cstart1;
                }
                if (s2) {
                    (*cl).cstart2 = ( (*tempsv).pos.start2 < (*cl).cstart2 ) ? (*tempsv).pos.start2 : (*cl).cstart2;
                }
                if (e1) {
                    (*cl).cend1 = ( (*tempsv).pos.end1 > (*cl).cend1 ) ? (*tempsv).pos.end1 : (*cl).cend1;
                }
                if (e2) {
                    (*cl).cend2 = ( (*tempsv).pos.end2 > (*cl).cend2 ) ? (*tempsv).pos.end2 : (*cl).cend2;
                }
            }
        }

        tol = ((*cl).tol * (numOfSVs + 1)) - (*sv).tolerance;
        (*cl).tol = tol / numOfSVs;
    }
    else {
        return -1;
    }

    return 1;
}

//Updates the number of edges of a cluster and degrees of its SVs according to the
//newly added or removed sv. Returns 1 if update is successful, 0 otherwise.
int updateCluster( Cluster *cl, SV *sv)
{
    SV *svptr, *crit, *rcrit;
    SV *tempsv;

    if ((*sv).clusterID == (*cl).id) //Update after addition
    {
        (*sv).clusterIndex = (*cl).svs.size() - 1;
        (*sv).degree = 0;
        for (int d = 0; d < (*cl).svs.size()-1; d++)
        {
            svptr = (*cl).svs[d];
	    if ((*sv).readID != (*svptr).readID)
            {

                if (assignEdge(sv, svptr, ovp) == 1)
                {
                    ((*cl).num_of_edges)++;
                    ((*sv).degree)++;
                    ((*svptr).degree)++;
                }
                else if ( (*svptr).degree < ceil( quasi_clique_lambda * (*cl).svs.size() ) )
                {
                    (*cl).crit.push_back(svptr);
                }
            }
        }

        if ( (*sv).degree < ceil( quasi_clique_lambda * (*cl).svs.size() ) )
            (*cl).crit.push_back(sv);

	if (updateBpAsCentroid(cl, sv, 0) != 1) {
	    return -1;
	}

	// If you plan to use this update, make sure that the overlap percentage and/or quasi clique parameters are tight!
/*	if (updateBpAsOutermost(cl, sv, 0) != 1) {
	    cout << "Start-end update failed." << endl;
	    return -1;
	}
*/
        //Update crit set
        for (int cr = 0; cr < (*cl).crit.size(); cr++)
        {
            crit = (*cl).crit[cr];
            if ( (*crit).degree >= ceil( quasi_clique_lambda * (*cl).svs.size() ) )
            {
                (*cl).crit.erase((*cl).crit.begin()+cr);
            }
        }

        return 1;
    }
    else
    {
        (*sv).clusterIndex = -1;
        (*sv).degree = 0;
        for (int d = 0; d < (*cl).svs.size(); d++)
        {
            svptr = (*cl).svs[d];
            if (assignEdge(sv, svptr, ovp))
            {
                ((*cl).num_of_edges)--;
                (*svptr).degree--;

                if ( ((*svptr).degree - 1) < ceil( quasi_clique_lambda * ((*cl).svs.size() - 2) ) )
                    (*cl).rCrit.push_back(svptr);
            }
        }

	if (updateBpAsCentroid(cl, sv, 1) != 1) {
	    return -1;
        }

        // If you plan to use this update, make sure that the overlap percentage and/or quasi clique parameters are tight!
/*      if (updateBpAsOutermost(cl, sv, 1) != 1) {
            cout << "Start-end update failed." << endl;
	    return -1;
        }
*/
        //Update rCrit set
        for (int rcr = 0; rcr < (*cl).rCrit.size(); rcr++)
        {
            rcrit = (*cl).rCrit[rcr];
            if ( ((*rcrit).degree - 1) >= ceil( quasi_clique_lambda * ((*cl).svs.size() - 2) ) )
            {
                (*cl).rCrit.erase((*cl).rCrit.begin()+rcr);
            }
        }

        return 1;
    }
}

//First, checks whether the chromosome of the new SV is the same as the cluster's.
//Then, checks whether the new SV has an edge to all the nodes in the critical
//addition set of the cluster. If it does, adds the new SV and returns 1. If not,
//returns 0. Returns -1 in case of errors.
int addSV( Cluster *cl, SV *newSV )
{
    if ( ((*cl).svs).empty() )
    {
        //Add pointer to the new SV to the cluster
        (*newSV).clusterID = (*cl).id;
        (*cl).svs.push_back(newSV);
        (*cl).chromosome = (*newSV).chromosome;

        //Update the cluster
        if (updateCluster(cl, newSV) == 1)
            return 1;
        else
        {
            (*cl).svs.pop_back();
            (*newSV).clusterID = -1;
            (*cl).chromosome = 0;
        }

        return 0;
    }

    if ((*cl).chromosome != (*newSV).chromosome)
    {
        return -1;
    }

    //The number of edges that the new SV contributes to the cluster
    int newEdges = 0;

    for (int i = 0; i < (*cl).crit.size(); i++)
    {
        if (assignEdge(newSV, (*cl).crit[i], ovp) == 0)
            return 0;
    }

    for (int i = 0; i < (*cl).svs.size(); i++)
    {
        if (assignEdge(newSV, (*cl).svs[i], ovp) == 1)
            newEdges++;
    }

    if ( ((*cl).num_of_edges + newEdges) < ceil( quasi_clique_gamma * ((*cl).svs.size() + 1) * (*cl).svs.size() / 2 ) )
        return 0;

    if ( newEdges < ceil( quasi_clique_lambda * (*cl).svs.size() ) )
        return 0;

    //Add pointer to the new SV to the cluster
    (*newSV).clusterID = (*cl).id;
    (*cl).svs.push_back(newSV);

    //Update the cluster
    if (updateCluster(cl, newSV) == 1)
        return 1;
    else
    {
        (*cl).svs.pop_back();
        (*newSV).clusterID = -1;
    }

    return 0;
}

//Removes an SV from a cluster.
int removeSV( Cluster *cl, SV *sv)
{
    SV *svptr;

    if ((*cl).svs.size() == 0 || (*sv).clusterID != (*cl).id)
        return -1;

    for (int rcr = 0; rcr < (*cl).rCrit.size(); rcr++)
    {
        if (assignEdge(sv, (*cl).rCrit[rcr], ovp) == 1)
            return 0;
    }

    if ( (*sv).degree > ( (*cl).num_of_edges - ceil( quasi_clique_gamma * ((*cl).svs.size() - 1) * ((*cl).svs.size() - 2) / 2 ) ) )
    {
        return 0;
    }

    //Remove SV
    (*sv).clusterID = -1;
    for (int d = 0; d < (*cl).svs.size(); d++)
    {
        svptr = (*cl).svs[d];
        if ( ((*svptr).pos.start1 == (*sv).pos.start1) && ((*svptr).pos.end1 == (*sv).pos.end1) )
        {
            (*cl).svs.erase((*cl).svs.begin()+d);
            break;
        }
    }

    //Update the cluster
    if (updateCluster(cl, sv) == 1)
        return 1;
    else
    {
        (*sv).clusterID = (*cl).id;
        (*cl).svs.push_back(sv);

        return 0;
    }
}

/*
 * Clusters the PacBio SVs using a quasi clique approximation algorithm.
 */
int clusterSVs( vector<SV> *svs, vector<Cluster> *clusters )
{
    if( (*svs).empty() )
    {
        return 0;
    }

    Cluster *cl;
    SV *sv;
    string svtype;

    Params *params = Params::getInstance();

    if ((*svs)[0].type == deletion) {
        quasi_clique_lambda = DEL_LAMBDA;
        quasi_clique_gamma = DEL_GAMMA;
	ovp = DEL_CLUSTER_OVERLAP_PERCENTAGE;
	minClusterSize = params->min_del_cluster_size;
	svtype = "Deletion";
    }
    else if ((*svs)[0].type == inversion) {
        quasi_clique_lambda = INV_LAMBDA;
        quasi_clique_gamma = INV_GAMMA;
	ovp = INV_CLUSTER_OVERLAP_PERCENTAGE;
	minClusterSize = params->min_inv_cluster_size;
	svtype = "Inversion";
    }
    else {
	quasi_clique_lambda = 0.5;
	quasi_clique_gamma = 0.6;
	minClusterSize = 1;
	ovp = 0.5;
	svtype = "none";
    }

    addCluster(clusters);

    cl = &((*clusters)[0]);
    sv = &((*svs)[0]);

    addSV(cl, sv);

    int c, fit;
    for (int d = 1; d < (*svs).size(); d++)
    {
        sv = &((*svs)[d]);

        fit = 0;
        for (c = 0; c < (*clusters).size(); c++)
        {
            cl = &((*clusters)[c]);

            if (addSV(cl, sv) == 1)
            {
                fit = 1;
                break;
            }
        }

        if (!fit)
        {
            addCluster(clusters);
            cl = &((*clusters)[c]);

            addSV(cl, sv);
        }
    }

    Cluster *cl_2;
    int repeat = 3;

    while (repeat != 0)
    {
        for (c = 0; c < (*clusters).size(); c++)
        {
            cl = &((*clusters)[c]);
            if ((*cl).svs.size() < minClusterSize )
            {
                for (int cc = 0; cc < (*clusters).size(); cc++)
                {
                    cl_2 = &((*clusters)[cc]);
                    if ((*cl).id != (*cl_2).id)
                    {
                        for (int d = 0; d < (*cl).svs.size(); d++)
                        {
                            sv = (*cl).svs[d];
                            if (addSV(cl_2, sv) == 1)
                                removeSV(cl, sv);
                        }
                    }
                }
            }
        }

        repeat--;
    }

    return 1;
}

void filterDups( vector<SV> *svsTemp, vector<int> sizeArr, vector<int> supportArr )
{
    for (int i = 0; i < (*svsTemp).size(); i++)
    {
	for (int j = i+1; j < (*svsTemp).size(); j++)
	{
	    if ((*svsTemp)[j].clusterID == -1) { continue; }

	    if (assignEdge(&((*svsTemp)[i]), &((*svsTemp)[j]), 0.5))
	    {
		if ((sizeArr[i]+supportArr[i]) < (sizeArr[j]+supportArr[j]))
		{
		    (*svsTemp)[i].clusterID = -1;
		    break;
		}
		else if ((sizeArr[j]+supportArr[j]) < (sizeArr[i]+supportArr[i]))
		{
		    (*svsTemp)[j].clusterID = -1;
		}
	    }
	}
    }
}

int filterSVs( vector<SV> *svs, vector<Cluster> *clusters, vector<SV> *finalSVs, const vector<short> *dV, const double meanReadDepth, const int currentChr, ofstream *out)
{
    SV clsv;
    vector<SV> svsTemp;
    vector<int> sizeArr, supportArr;

    int found, o = -2, con = -1, valids = 0;
    uint64_t start1, start2, end1, end2;
    int lowclsize, lowmatchcount;
    string svtype;

    double readDepthPerc, meanSize=0, meanSupport=0;
    double stdSize=0, stdSupport=0;
    int totalSV=0;

    Params *params = Params::getInstance();

    if ((*clusters).empty() || (*svs).empty()) {
	return 0;
    }

    if ((*svs)[0].type == deletion)
    {
	minClusterSize = params->min_del_cluster_size;
        minSupport = params->min_del_support;
	depthThreshold = DEL_DEPTH_THRESHOLD;
	ovp = DEL_CLUSTER_OVERLAP_PERCENTAGE;
	svtype = "Deletion";
    }
    else if ((*svs)[0].type == inversion)
    {
	minClusterSize = params->min_inv_cluster_size;
        minSupport = params->min_inv_support;
	depthThreshold = INV_DEPTH_THRESHOLD;
	ovp = INV_CLUSTER_OVERLAP_PERCENTAGE;
	svtype ="Inversion";
    }
    else{
        minClusterSize = 1;
        minSupport = 1;
	depthThreshold = 3;
	ovp = 0.5;
	svtype = "none";
    }

    lowclsize = 0;
    lowmatchcount = 0;

    clsv.type = (*svs)[0].type;
    for (int c = 0; c < (*clusters).size(); c++)
    {
        found = 0;
	clsv.chromosome = (*clusters)[c].chromosome;
	clsv.pos.start1 = (*clusters)[c].cstart1;
	clsv.pos.start2 = (*clusters)[c].cstart2;
	clsv.pos.end1 = (*clusters)[c].cend1;
	clsv.pos.end2 = (*clusters)[c].cend2;
	clsv.tolerance = (*clusters)[c].tol;
	clsv.clusterID = (*clusters)[c].id;

        for (int d = 0; d < (*svs).size(); d++)
        {
	    o = assignEdge( &clsv, &((*svs)[d]), ovp );
	    if (o == 1)
	    {
		found++;
	    }
	}

	if (found == 0) { continue; }

	if ( ( (*clusters)[c].svs.size() < minClusterSize && found < minSupport ) && ( ((*clusters)[c].svs.size() + found) < ((minClusterSize + minSupport) / 2) ) )
	{
	    (*clusters)[c].id = -1;

	    lowclsize++;
	    lowmatchcount++;

	    continue;
	}

        readDepthPerc = 0;
        for( int i = clsv.pos.start1; i < clsv.pos.end2; i++ )
            readDepthPerc += (*dV)[i];

	readDepthPerc = readDepthPerc / ( meanReadDepth * (double)(clsv.pos.end2-clsv.pos.start1) );

	if (readDepthPerc <= depthThreshold) {
	    svsTemp.push_back(clsv);
	    sizeArr.push_back((*clusters)[c].svs.size());
	    supportArr.push_back(found);
	}
    }

    filterDups(&svsTemp, sizeArr, supportArr);

    for( int i = 0; i < svsTemp.size(); i++)
    {
	if (svsTemp[i].clusterID == -1) {
            continue; }

	totalSV++;
	meanSize += sizeArr[i];
	meanSupport += supportArr[i];
    }

    meanSize = meanSize/totalSV;
    meanSupport = meanSupport/totalSV;

    for( int i = 0; i < svsTemp.size(); i++)
    {
        if (svsTemp[i].clusterID == -1) {
            continue; }

        stdSize += (meanSize - sizeArr[i])*(meanSize - sizeArr[i]);
        stdSupport += (meanSupport - supportArr[i])*(meanSupport - supportArr[i]);
    }

    stdSize = sqrt(stdSize);
    stdSupport = sqrt(stdSupport);

    for( int i = 0; i < svsTemp.size(); i++)
    {
	if (svsTemp[i].clusterID == -1) {
	    continue; }

	if (svsTemp[i].type == deletion && (sizeArr[i] < meanSize) && (supportArr[i] < meanSupport)) {
	    lowclsize++;
	    lowmatchcount++;
	    continue;
	}
	else if ( svsTemp[i].type == inversion && ( (sizeArr[i] < meanSize) || (supportArr[i] < meanSupport) ) ) {
	    if (sizeArr[i] < meanSize) { lowclsize++; }
	    if (supportArr[i] < meanSupport) { lowmatchcount++; }
	    continue;
	}

        valids++;

        //Store the final SV
        (*finalSVs).push_back(svsTemp[i]);

        if ((*out).is_open())
        {
            if (clsv.type == deletion){
                (*out) << clsv.chromosome << '\t' << clsv.pos.start1 << '\t' << clsv.pos.end1 << '\n';
            }
            else {
                (*out) << clsv.chromosome << '\t' << clsv.pos.start1 << '\t' << clsv.pos.start2 << '\t' << clsv.pos.end1
                       << '\t' << clsv.pos.end2 << '\n';
            }
        }
        else { return 0; }
    }

    return 1;
}
