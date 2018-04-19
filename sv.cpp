#include "sv.h"

Params* Params::instance = NULL;

int validateRead( Read *rd, int currentChr, sonic *snc, int type )
{
    int end;
    int inDupRegion, inSatRegion;
    int push;
    int filterSegDups, filterSats;

    end = (*rd).pos + referenceSpan(rd);

    inDupRegion = sonic_is_segmental_duplication(snc, snc->chromosome_names[currentChr-1], (*rd).pos, end);
    inSatRegion = sonic_is_satellite(snc, snc->chromosome_names[currentChr-1], (*rd).pos, end);

    push = 0;

    if ( type == 0 ) {
        filterSegDups = DEL_FILTER_SEGMENTAL_DUPLICATIONS;
        filterSats = DEL_FILTER_SATELLITES;
    }
    else {
        filterSegDups = INV_FILTER_SEGMENTAL_DUPLICATIONS;
        filterSats = INV_FILTER_SATELLITES;
    }

    if ( !filterSegDups && !filterSats )
    {
        push = 1;
    }
    else if ( !filterSegDups && filterSats )
    {
        if (inSatRegion == 0) { push = 1; }
    }
    else if ( filterSegDups && !filterSats )
    {
        if (inDupRegion == 0) { push = 1; }
    }
    else
    {
        if (inSatRegion == 0 && inDupRegion == 0) { push = 1; }
    }

    return push;
}

int percOvlp( const int sv1start, const int sv2start, const int sv1end, const int sv2end, double overlapPcg )
{
    int nonOverlap = 0;
    if ( sv1start > sv2start )
    {
        nonOverlap += sv1start - sv2start;
    }
    else
    {
        nonOverlap += sv2start - sv1start;
    }

    if ( sv1end > sv2end )
    {
        nonOverlap += sv1end - sv2end;
    }
    else
    {
        nonOverlap += sv2end - sv1end;
    }

    int sv1_len = sv1end - sv1start;
    int sv2_len = sv2end - sv2start;
    int total = sv1_len + sv2_len;

    double percentageOverlap = ((double)(total - nonOverlap) / 2.0) / (double)sv1_len;
    if ( percentageOverlap > overlapPcg )
    {
        percentageOverlap = ((double)(total - nonOverlap) / 2.0) / (double)sv2_len;

        if ( percentageOverlap > overlapPcg )
        {
            return 1;
        }
    }

    return 0;
}

int toleranceOvlp( const SV *s1, const SV *s2 )
{
    if ((*s1).chromosome != (*s2).chromosome)
    {
        return 0;
    }

    //Tolerance interval for when comparing two SVs is determined as the average of
    //their individual tolerance intervals.
    int toleranceInterval = ((*s1).tolerance + (*s2).tolerance) / 2;

    if ( ((*s1).pos.start1 <= (*s2).pos.start1 + toleranceInterval)
         && ((*s1).pos.start1 >= (*s2).pos.start1 - toleranceInterval)
         && ((*s1).pos.end1 <= (*s2).pos.end1 + toleranceInterval)
         && ((*s1).pos.end1 >= (*s2).pos.end1 - toleranceInterval) )
        return 1;

    return 0;
}

int delOvlp( const SV *s1, const SV *s2, double overlapPcg )
{
    return percOvlp((*s1).pos.start1, (*s2).pos.start1, (*s1).pos.end1, (*s2).pos.end1, overlapPcg);

}

int invOvlp( const SV *s1, const SV *s2, double overlapPcg )
{
//    return ( percOvlp((*s1).pos.start1, (*s2).pos.start1, (*s1).pos.end1, (*s2).pos.end1, overlapPcg)
//             & percOvlp((*s1).pos.start2, (*s2).pos.start2, (*s1).pos.end2, (*s2).pos.end2, overlapPcg) );

    return ( percOvlp( ((*s1).pos.start1+(*s1).pos.start2)/2, ((*s2).pos.start1+(*s2).pos.start2)/2,
	   ((*s1).pos.end1+(*s1).pos.end2)/2 ,((*s2).pos.end1+(*s2).pos.end2)/2, overlapPcg ) );
}

int assignEdge( const SV *s1, const SV *s2, double ovp )
{
    if ((*s1).chromosome != (*s2).chromosome) {
        return 0;
    }

    if ((*s1).type == deletion) {
        return delOvlp(s1, s2, ovp);
    }

    if ((*s1).type == inversion) {
        return invOvlp(s1, s2, ovp);
    }

    else { return -1; }
}

int mergeInv( const SV *s1, const SV *s2)
{
    return toleranceOvlp(s1, s2);
}

int delCigar( const Read *rd, uint64_t *start, uint64_t *len, uint64_t *tol )
{
    int op, op_len, interference, m_intf, total_deletion_length, processing;
    uint64_t refStart;

    vector<uint64_t> del_start;
    vector<int> del_len;
    int read_length = readLength(rd);

    Params *params = Params::getInstance();

    *tol = 0.15 * read_length;

    refStart = 0;
    interference = 0;
    m_intf = 0;
    total_deletion_length = 0;

    processing = 0;
    for (int i = 0; i < (*rd).ncigar; i++)
    {
        op = bam_cigar_op((*rd).cigar[i]);

        if (op == 2)
        {
//	    cout << "deletion" << endl;

            op_len = bam_cigar_oplen((*rd).cigar[i]);
            if( !processing )
            {
//		cout << "processing = 0, ";
                if( op_len < MIN_DEL_STARTING_LEN )
                    continue;
                *start = (*rd).pos + refStart;
//		cout << "start: " << *start << endl;
                processing = 1;
            }

            total_deletion_length += m_intf;
            m_intf = 0;
            interference = 0;

            total_deletion_length += op_len;
            *len = total_deletion_length;

            refStart += op_len;

//	    cout << "deletion len: " << op_len << ", totalDL: " << total_deletion_length << ", refStart: " << refStart << endl;
        }

        if (op == 0)
        {
            op_len = bam_cigar_oplen((*rd).cigar[i]);
            refStart += op_len;

//	    cout << "match. match len: " << op_len << ", refStart: " << refStart << endl;

            if( processing )
            {
                interference += op_len;
                m_intf += op_len;
//		cout << "processing = 1, interference: " << interference;
                if( interference > INTERFERING_BPS )
                {
//		    cout << ", totalDL: " << total_deletion_length << endl;
                    if( (total_deletion_length >= params->min_del_len) && (total_deletion_length <= MAX_DEL_LEN) )
                        return 1;

                    interference = 0;
                    m_intf = 0;
                    total_deletion_length = 0;
                    processing = 0;
                }
//		cout << endl;
            }
        }

        if (op == 1)
        {
//	    cout << "insertion" << endl;
            if( processing )
            {
                op_len = bam_cigar_oplen((*rd).cigar[i]);
                interference += op_len;
//		cout << "processing = 1, insertion len: " << op_len << ", interference: " << interference;
                if( interference > INTERFERING_BPS )
                {
//		    cout << ", totalDL: " << total_deletion_length << endl;
                    if( (total_deletion_length >= params->min_del_len) && (total_deletion_length <= MAX_DEL_LEN) )
                        return 1;

                    interference = 0;
                    m_intf = 0;
                    total_deletion_length = 0;
                    processing = 0;
                }
//		cout << endl;
            }
        }
    }

//  cout << ", totalDL: " << total_deletion_length << endl;
    if( (total_deletion_length >= params->min_del_len) && (total_deletion_length <= MAX_DEL_LEN) )
        return 1;
    else
        return 0;
}

void detectLowDepth( vector<SV> *all_deletions, const vector<SV> *il_deletions, const int currentChr, sonic *snc, const vector<short> *depthVector, ofstream *out )
{
    SV del, sv;
    int del_len;
    int inGapRegion, inDupRegion, inSatRegion;
    int push, found;
    uint64_t start1, start2, end1, end2;

    Params *params = Params::getInstance();

    for (unsigned i = 0; i < (*depthVector).size(); i++)
    {
        if ((*depthVector)[i] <= 1)
        {
            del.pos.start1 = i;
            del.pos.start2 = i;

            while ((i < (*depthVector).size()) && ((*depthVector)[i] <= 1))
            {
                i++;
            }

            del.pos.end1 = i-1;
            del.pos.end2 = i-1;

	    inDupRegion = sonic_is_segmental_duplication(snc, snc->chromosome_names[currentChr-1], del.pos.start1, del.pos.end2);
	    inSatRegion = sonic_is_satellite(snc, snc->chromosome_names[currentChr-1], del.pos.start1, del.pos.end2);
	    inGapRegion = sonic_is_gap(snc, snc->chromosome_names[currentChr-1], del.pos.start1, del.pos.end2);

	    if (inGapRegion) {
		continue;
	    }

	    push = 0;
	    if ( !INV_FILTER_SEGMENTAL_DUPLICATIONS && !INV_FILTER_SATELLITES ) {
		push = 1; }
	    else if ( !INV_FILTER_SEGMENTAL_DUPLICATIONS && INV_FILTER_SATELLITES ) {
		if (inSatRegion == 0) { push = 1; } }
	    else if ( INV_FILTER_SEGMENTAL_DUPLICATIONS && !INV_FILTER_SATELLITES ) {
		if (inDupRegion == 0) { push = 1; } }
	    else {
		if (inSatRegion == 0 && inDupRegion == 0) { push = 1; } }

	    if (!push) { continue; }

            del_len = del.pos.end1 - del.pos.start1;
	    del.chromosome = currentChr;
	    del.tolerance = del_len*0.15;
	    del.type = deletion;
	    del.readID = "noread";
            if (del_len >= params->min_del_len && del_len <= MAX_DEL_LEN)
            {
		found = 0;
	        for (int d = 0; d < (*il_deletions).size(); d++)
	        {
	            int o = assignEdge( &del, &((*il_deletions)[d]), DEL_VERIFICATION_OVERLAP_PERCENTAGE );
	            if (o == 1)
	            {
	                found++;

			sv = (*il_deletions)[d];
			if( sv.pos.start1 > del.pos.start1 ) {
			    del.pos.start1 = sv.pos.start1;
			}
			if( sv.pos.start2 > del.pos.start2 ) {
			    del.pos.start2 = sv.pos.start2;
			}
			if( sv.pos.end1 < del.pos.end1 ) {
			   del.pos.end1 = sv.pos.end1;
			}
			if( sv.pos.end2 < del.pos.end2 ) {
			    del.pos.end2 = sv.pos.end2;
			}

	            }
	        }

//		if (found < DEL_ILLUMINA_MATCH_COUNT)
		if (found < 5)
	        {
	            continue;
	        }

                (*all_deletions).push_back(del);
            }
        }
    }
}

int verifyInversions( const vector<SV> *forwardSplits, const vector<SV> *reverseSplits, vector<SV> *inversions, double ovpc, int merge )
{
    SV inv_f, inv_r, newInv;

    cout << "Verifying inversions.\nNumber of forward splits: " << (*forwardSplits).size() << endl;
    cout << "Number of reverse splits: " << (*reverseSplits).size() << endl;

    if (merge == 0)
    {
	for (unsigned i = 0; i < (*forwardSplits).size(); i++)
	{
	    (*inversions).push_back((*forwardSplits)[i]);
	}

	for (unsigned i = 0; i < (*reverseSplits).size(); i++)
        {
            (*inversions).push_back((*reverseSplits)[i]);
        }

	return 1;
    }

    for (unsigned i = 0; i < (*forwardSplits).size(); i++)
    {
        for (unsigned j = 1; j < (*reverseSplits).size(); j++)
        {
            inv_f = (*forwardSplits)[i];
            inv_r = (*reverseSplits)[j];

//          if( mergeInv(&inv_f, &inv_r) == 1 )
	    if( delOvlp(&inv_f, &inv_r, ovpc) == 1)
            {
                // add inversion to inversions vector
                newInv.chromosome = inv_f.chromosome;
                newInv.tolerance = (inv_f.tolerance + inv_r.tolerance) / 2;
		newInv.type = inversion;
                newInv.readID = inv_f.readID + "&" + inv_r.readID;

                if (inv_f.pos.start1 < inv_r.pos.start1)
                {
                    newInv.pos.start1 = inv_f.pos.start1;
                    newInv.pos.start2 = inv_r.pos.start1;
                }
                else
                {
                    newInv.pos.start1 = inv_r.pos.start1;
                    newInv.pos.start2 = inv_f.pos.start1;
                }

                if (inv_f.pos.end1 < inv_r.pos.end1)
                {
                    newInv.pos.end1 = inv_f.pos.end1;
                    newInv.pos.end2 = inv_r.pos.end1;
                }
                else
                {
                    newInv.pos.end1 = inv_r.pos.end1;
                    newInv.pos.end2 = inv_f.pos.end1;
                }

//		cout << "Pushing back inversion..." << endl;
                (*inversions).push_back(newInv);
            }
        }
    }

    return 1;
}

int findIlluminaDel( vector<SV> *deletions, const bam1_t *aln, const int currentChr, sonic *snc )
{
    int deletionLength;
    int inGapRegion, inSatRegion, inDupRegion;

    SV del, altDel;

    Params *params = Params::getInstance();

    del.chromosome = aln->core.tid+1;
    del.pos.start1 = bam_endpos(aln);
    del.pos.start2 = bam_endpos(aln);
    del.pos.end1 = aln->core.mpos - 1;
    del.pos.end2 = aln->core.mpos - 1;
    deletionLength = del.pos.end1 - del.pos.start1 + 1;

    if ( (deletionLength >= params->min_del_len) && (deletionLength <= MAX_DEL_LEN) )
    {
	inDupRegion = sonic_is_segmental_duplication(snc, snc->chromosome_names[currentChr-1], del.pos.start1, del.pos.end2);
        inSatRegion = sonic_is_satellite(snc, snc->chromosome_names[currentChr-1], del.pos.start1, del.pos.end2);
        inGapRegion = sonic_is_gap(snc, snc->chromosome_names[currentChr-1], del.pos.start1, del.pos.end1);

        if (inGapRegion) { return 0; }

    int push = 0;

    if ( !DEL_FILTER_SEGMENTAL_DUPLICATIONS && !DEL_FILTER_SATELLITES )
    {
        push = 1;
    }
    else if ( !DEL_FILTER_SEGMENTAL_DUPLICATIONS && DEL_FILTER_SATELLITES )
    {
        if (inSatRegion == 0) { push = 1; }
    }
    else if ( DEL_FILTER_SEGMENTAL_DUPLICATIONS && !DEL_FILTER_SATELLITES )
    {
        if (inDupRegion == 0) { push = 1; }
    }
    else
    {
        if (inSatRegion == 0 && inDupRegion == 0) { push = 1; }
    }
        del.readID = (std::string) bam_get_qname(aln);
        del.type = deletion;

	if (push) {
            (*deletions).push_back(del);
        }
    }
    return 1;
}

int findIlluminaInv( vector<SV> *inversions, vector<SV> *forwardPairs, vector<SV> *reversePairs, const bam1_t *aln, const int currentChr, sonic *snc )
{
    int inversionLength;
    int inGapRegion, inSatRegion, inDupRegion;
    int inReverseStrand, mateInReverseStrand;

    SV inv, altInv;

    Params *params = Params::getInstance();

    inv.chromosome = aln->core.tid+1;
    inv.pos.start1 = bam_endpos(aln) + 1;
    inv.pos.start2 = bam_endpos(aln) + 1;
    inv.pos.end1 = aln->core.mpos + mateReferenceSpan(aln);
    inv.pos.end2 = aln->core.mpos + mateReferenceSpan(aln);
    inversionLength = inv.pos.end2 - inv.pos.start1 + 1;

    if ( (inversionLength >= params->min_inv_len) && (inversionLength <= MAX_INV_LEN) )
    {
        inGapRegion = sonic_is_gap(snc, snc->chromosome_names[currentChr-1], inv.pos.start1, inv.pos.end2);

        if (inGapRegion) {
	    return 0;
	}

	inDupRegion = sonic_is_segmental_duplication(snc, snc->chromosome_names[currentChr-1], inv.pos.start1, inv.pos.end2);
            inSatRegion = sonic_is_satellite(snc, snc->chromosome_names[currentChr-1], inv.pos.start1, inv.pos.end2);

            int push = 0;
            if ( !INV_FILTER_SEGMENTAL_DUPLICATIONS && !INV_FILTER_SATELLITES ) {
                push = 1; }
            else if ( !INV_FILTER_SEGMENTAL_DUPLICATIONS && INV_FILTER_SATELLITES ) {
                if (inSatRegion == 0) { push = 1; } }
            else if ( INV_FILTER_SEGMENTAL_DUPLICATIONS && !INV_FILTER_SATELLITES ) {
                if (inDupRegion == 0) { push = 1; } }
            else {
                if (inSatRegion == 0 && inDupRegion == 0) { push = 1; } }

            if (!push) { return 0; }

        inv.readID = (std::string) bam_get_qname(aln);
	inv.type = inversion;

        inReverseStrand = aln->core.flag&BAM_FREVERSE;
        mateInReverseStrand = aln->core.flag&BAM_FMREVERSE;

	if ( (!inReverseStrand && !mateInReverseStrand) || (inReverseStrand && mateInReverseStrand) )
	{
	     (*inversions).push_back(inv);
	}
    }

    return 1;
}

int findPacbioDel( vector<SV> *deletions, const Read *cur_rd, const Read *next_rd, const int currentChr, sonic *snc )
{
    SV del;

    vector<int> currentSClengths, nextSClengths;

    int currentReadStart, currentReadEnd, nextReadStart, nextReadEnd;	//Read variables
    uint64_t currentReferenceStart, currentReferenceEnd, nextReferenceStart, nextReferenceEnd, toleranceInterval;  //Reference variables
    int currentReadLength, nextReadLength;
    int inGapRegion, inSatRegion, inDupRegion;

    Params *params = Params::getInstance();

    softClipLengths(cur_rd, &currentSClengths);

    currentReadStart = currentSClengths[0];
    currentReadLength = readLength(cur_rd);
    currentReadEnd = currentReadStart + currentReadLength - 1;

    softClipLengths(next_rd, &nextSClengths);

    nextReadStart = nextSClengths[0];
    nextReadLength = readLength(next_rd);
    nextReadEnd = nextReadStart + nextReadLength - 1;

    toleranceInterval = ( currentReadLength + nextReadLength ) * 0.15;

    // Check split condition
    if ( ( (nextReadStart <= (currentReadEnd + toleranceInterval)) && (nextReadStart >= (currentReadEnd - toleranceInterval)) )
         || ( (currentReadStart <= (nextReadEnd + toleranceInterval)) && (currentReadStart >= (nextReadEnd - toleranceInterval)) ) )
//	if ( ( (nextReadStart <= (currentReadEnd + 2)) && (nextReadStart >= currentReadEnd) )
//		 || ( (currentReadStart <= (nextReadEnd + 2)) && (currentReadStart >= nextReadEnd) ) )
    {
        currentReferenceStart = (*cur_rd).pos;
        currentReferenceEnd = currentReferenceStart + referenceSpan(cur_rd) - 1;

        nextReferenceStart = (*next_rd).pos;        // inclusive
        nextReferenceEnd = nextReferenceStart + referenceSpan(next_rd) - 1;     // inclusive

        if ( ((nextReferenceStart - currentReferenceEnd - 1) >= params->min_del_len) && ((nextReferenceStart - currentReferenceEnd - 1) <= MAX_DEL_LEN) )
        {
            del.pos.start1 = currentReferenceEnd+1;
            del.pos.start2 = currentReferenceEnd+1;
            del.pos.end1 = nextReferenceStart-1;
            del.pos.end2 = nextReferenceStart-1;
        }
        else if ( ((currentReferenceStart - nextReferenceEnd - 1) >= params->min_del_len) && ((currentReferenceStart - nextReferenceEnd - 1) <= MAX_DEL_LEN) )
        {
            del.pos.start1 = nextReferenceEnd+1;
            del.pos.start2 = nextReferenceEnd+1;
            del.pos.end1 = currentReferenceStart-1;
            del.pos.end2 = currentReferenceStart-1;
        }
        else
        {
            return 0;
        }

        inGapRegion = sonic_is_gap(snc, snc->chromosome_names[currentChr-1], del.pos.start1, del.pos.end2);

        if (inGapRegion) { return 0; }

	inDupRegion = sonic_is_segmental_duplication(snc, snc->chromosome_names[currentChr-1], del.pos.start1, del.pos.end2);
        inSatRegion = sonic_is_satellite(snc, snc->chromosome_names[currentChr-1], del.pos.start1, del.pos.end2);

    int push = 0;

    if ( !DEL_FILTER_SEGMENTAL_DUPLICATIONS && !DEL_FILTER_SATELLITES )
    {
        push = 1;
    }
    else if ( !DEL_FILTER_SEGMENTAL_DUPLICATIONS && DEL_FILTER_SATELLITES )
    {
        if (inSatRegion == 0) { push = 1; }
    }
    else if ( DEL_FILTER_SEGMENTAL_DUPLICATIONS && !DEL_FILTER_SATELLITES )
    {
        if (inDupRegion == 0) { push = 1; }
    }
    else
    {
        if (inSatRegion == 0 && inDupRegion == 0) { push = 1; }
    }
	if(!push) { return 0; }

        del.chromosome = (*cur_rd).rname;
        del.tolerance = toleranceInterval;
        del.readID = (*cur_rd).qname;
        del.type = deletion;

//	cout << "Pushing back deletion: " << del.readID << endl;
        (*deletions).push_back(del);
    }

    return 1;
}

int findPacbioInv( vector<SV> *forwardSplits, vector<SV> *reverseSplits, const Read *cur_rd, const Read *next_rd, const int currentChr, sonic *snc )
{
    SV inv;

    vector<int> currentSClengths, nextSClengths;

    int currentReadStart, currentReadEnd, nextReadStart, nextReadEnd;	//Read variables
    uint64_t currentReferenceStart, currentReferenceEnd, nextReferenceStart, nextReferenceEnd, toleranceInterval;  //Reference variables
    int currentReadLength, nextReadLength;
    int inGapRegion, inSatRegion, inDupRegion;
    int opposite;

    Params *params = Params::getInstance();

    softClipLengths(cur_rd, &currentSClengths);

    currentReadStart = currentSClengths[0];
    currentReadLength = readLength(cur_rd);
    currentReadEnd = currentReadStart + currentReadLength - 1;

    softClipLengths(next_rd, &nextSClengths);

    nextReadStart = nextSClengths[0];
    nextReadLength = readLength(next_rd);
    nextReadEnd = nextReadStart + nextReadLength - 1;

    toleranceInterval = ( currentReadLength + nextReadLength ) * 0.15;

    // Check split condition
    if ( ( (nextReadStart <= (currentReadEnd + toleranceInterval)) && (nextReadStart >= (currentReadEnd - toleranceInterval)) )
         || ( (currentReadStart <= (nextReadEnd + toleranceInterval)) && (currentReadStart >= (nextReadEnd - toleranceInterval)) ) )
    {
        currentReferenceStart = (*cur_rd).pos;
        currentReferenceEnd = currentReferenceStart + referenceSpan(cur_rd) - 1;

        nextReferenceStart = (*next_rd).pos;        // inclusive
        nextReferenceEnd = nextReferenceStart + referenceSpan(next_rd) - 1;     // inclusive

        if ( ((nextReferenceEnd - currentReferenceEnd - 1) >= params->min_inv_len) && ((nextReferenceEnd - currentReferenceEnd - 1) <= MAX_INV_LEN) )
        {
            inv.pos.start1 = currentReferenceEnd+1;
	    inv.pos.start2 = currentReferenceEnd+1;
            inv.pos.end1 = nextReferenceEnd-1;
	    inv.pos.end2 = nextReferenceEnd-1;

            opposite = 0;
        }
        else if ( ((currentReferenceEnd - nextReferenceEnd - 1) >= params->min_inv_len) && ((currentReferenceEnd - nextReferenceEnd - 1) <= MAX_INV_LEN) )
        {
            inv.pos.start1 = nextReferenceEnd+1;
	    inv.pos.start2 = nextReferenceEnd+1;
            inv.pos.end1 = currentReferenceEnd-1;
	    inv.pos.end2 = currentReferenceEnd-1;

            opposite = 1;
        }
        else
        {
            return -1;
        }

        inGapRegion = sonic_is_gap(snc, snc->chromosome_names[currentChr-1], inv.pos.start1, inv.pos.end2);

        if (inGapRegion) {
	    return 0;
	}

	inDupRegion = sonic_is_segmental_duplication(snc, snc->chromosome_names[currentChr-1], inv.pos.start1, inv.pos.end2);
            inSatRegion = sonic_is_satellite(snc, snc->chromosome_names[currentChr-1], inv.pos.start1, inv.pos.end2);

            int push = 0;
            if ( !INV_FILTER_SEGMENTAL_DUPLICATIONS && !INV_FILTER_SATELLITES ) {
                push = 1; }
            else if ( !INV_FILTER_SEGMENTAL_DUPLICATIONS && INV_FILTER_SATELLITES ) {
                if (inSatRegion == 0) { push = 1; } }
            else if ( INV_FILTER_SEGMENTAL_DUPLICATIONS && !INV_FILTER_SATELLITES ) {
                if (inDupRegion == 0) { push = 1; } }
            else {
                if (inSatRegion == 0 && inDupRegion == 0) { push = 1; } }

            if (!push) { return 0; }

        inv.chromosome = (*cur_rd).rname;
        inv.tolerance = toleranceInterval;
	inv.type = inversion;
        inv.readID = (*cur_rd).qname;

        if ( ( ( ((*cur_rd).flag) & BAM_FREVERSE ) == 0 ) && ( ( ((*next_rd).flag) & BAM_FREVERSE ) != 0 ) )
        {
            if( !opposite ){
                (*forwardSplits).push_back(inv);
            }
            else{ (*reverseSplits).push_back(inv); }
        }
        else if ( ( ( ((*cur_rd).flag) & BAM_FREVERSE ) != 0 ) && ( ( ((*next_rd).flag) & BAM_FREVERSE ) == 0 ) )
        {
            if( !opposite ){
                (*reverseSplits).push_back(inv);
            }
            else{ (*forwardSplits).push_back(inv); }
        }
    }

    currentSClengths.clear();
    nextSClengths.clear();

    return 1;
}


int findPacbioSVs( Unmmap *candidates, vector<SV> *deletions, vector<SV> *inversions, const int currentChr, sonic *snc )
{
    if( (*candidates).empty() )
    {
        return -1;
    }

    vector<SV> forwardSplits, reverseSplits;

    int curInReverseStrand, nextInReverseStrand;	//Flag variables
    int curBucketPos = 0;

    string cur_key, next_key;
    Read *cur_rd, *next_rd;

    for (auto currentRead = (*candidates).begin(); currentRead != (*candidates).end(); currentRead++)
    {
        cur_key = currentRead->first;

        if (curBucketPos == (*candidates).bucket_size((*candidates).bucket(cur_key)) - 1)
        {
            curBucketPos = 0;
            continue;
        }

        cur_rd = &(currentRead->second);
        if( currentChr != (*cur_rd).rname )
        {
            cout << (*cur_rd).rname;
            continue;
        }

        curInReverseStrand = (*cur_rd).flag&BAM_FREVERSE;
        if (curInReverseStrand > 0) {
            curInReverseStrand = 1;
	}

        auto nextRead = (*candidates).begin((*candidates).bucket(cur_key));
        for( int i = 0; i < curBucketPos+1; i++ )
        {
            nextRead++;
        }

        for ( ; nextRead != (*candidates).end((*candidates).bucket(cur_key)); ++nextRead)
        {
            next_key = nextRead->first;
            if (cur_key.compare(next_key) == 0)
            {
                next_rd = &(nextRead->second);

                nextInReverseStrand = (*next_rd).flag&BAM_FREVERSE;
                if (nextInReverseStrand > 0) {
                    nextInReverseStrand = 1;
		}

                if ( (curInReverseStrand == 1 && nextInReverseStrand == 1) || (curInReverseStrand != 1 && nextInReverseStrand != 1) )
                {
                    findPacbioDel (deletions, cur_rd, next_rd, currentChr, snc);
                }

                if ( (curInReverseStrand == 1 && nextInReverseStrand != 1) || (curInReverseStrand != 1 && nextInReverseStrand == 1) )
                {
                     findPacbioInv (&forwardSplits, &reverseSplits, cur_rd, next_rd, currentChr, snc);
                }
            }
        }

	curBucketPos++;
    }

    verifyInversions(&forwardSplits, &reverseSplits, inversions, DEL_VERIFICATION_OVERLAP_PERCENTAGE, 1);

    return 1;
}
