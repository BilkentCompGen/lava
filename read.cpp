//
// Created by Ezgi Ebren on 22.03.2018.
//

/*
 * Input: Read *rd - A read pointer
 *	  int numofargs - The number of cigar operations passed as arguments to the function
 *	  ... - Cigar operations that you want to count (int)
 *
 * Returns: A list (integer) of total counts of each cigar operation given in the arguments
 */

#include "read.h"

vector<int> totalOpCount( const Read *rd, int numofargs, ... )
{
    va_list cigarOps;

    vector<int> opCounts;
    int op, op_len;

    for( int i = 0; i < numofargs; i++ ){
        opCounts.push_back(0);
    }

    for( int j = 0; j < (*rd).ncigar; j++ )
    {
        op = bam_cigar_op((*rd).cigar[j]);

        va_start (cigarOps, numofargs);
        for( int i = 0; i < numofargs; i++ )
        {
            if( op == va_arg(cigarOps, int) )
            {
                op_len = bam_cigar_oplen((*rd).cigar[j]);
                if( op_len < 0 )
                {
                    cout << op_len << ", ";
                    op_len = op_len * (-1);
                }
                opCounts[i] += op_len;
            }
        }
        va_end (cigarOps);
    }

    return opCounts;
}

/*
 * Finds the length (in terms of number of bases) of the soft clipped (sc) regions at the beginning and end of a read.
 * Outputs a pair of integers, first one denoting the length of the sc at the beginning of the read and second one
 * denoting the length of the sc at the end of the read.
 */
void softClipLengths( const Read *rd, vector<int> *sclengths )
{
    int op_left, op_len_left, op_right, op_len_right;

    op_left = bam_cigar_op((*rd).cigar[0]);
    op_right = bam_cigar_op((*rd).cigar[(*rd).ncigar-1]);

    if (op_left == 4)
    {
        op_len_left = bam_cigar_oplen((*rd).cigar[0]);
        (*sclengths).push_back(op_len_left);
    }
    else
    {
        (*sclengths).push_back(0);
    }

    if (op_right == 4)
    {
        op_len_right = bam_cigar_oplen((*rd).cigar[(*rd).ncigar-1]);
        (*sclengths).push_back(op_len_right);
    }
    else
    {
        (*sclengths).push_back(0);
    }
}

/*
 * Computes read reference span (#M + #D)
 */
int referenceSpan( const Read *rd )
{
    vector<int> ops;
    ops = totalOpCount(rd, 2, BAM_CMATCH, BAM_CDEL);

    return (ops[0] + ops[1]);
}

/*
 * Computes read length (#M + #I)
 */
int readLength( const Read *rd )
{
    vector<int> ops;
    ops = totalOpCount(rd, 2, BAM_CMATCH, BAM_CINS);

    return (ops[0] + ops[1]);
}

/*
 * Computes illumina mate reference span
 */
int mateReferenceSpan( const bam1_t *aln )
{
    int referenceSpan = 0;
    char *mateCigar;
    string op_len = "";

    uint8_t *mc = bam_aux_get(aln, "MC");
    if (mc)
    {
        mateCigar = bam_aux2Z( mc );
        for (char *c = mateCigar; *c; ++c)
        {
            if (int(*c) >= 48 && int(*c) <= 57)
            {
                op_len += *c;
            }

            if (*c == 'S' || *c == 'I')
                op_len = "";

            if (*c == 'M' || *c == 'D')
            {
                referenceSpan += stoi(op_len);
                op_len = "";
            }
        }
    }
    else
        referenceSpan = 0;

    return referenceSpan;
}
/*
 * Computes alternative reference span (#M + #D)
 */
int altRefSpan(const char *altCigar)
{
    int referenceSpan = 0;
    string op_len = "";

    for (char *c = (char*)altCigar; *c; ++c)
    {
        if (int(*c) >= 48 && int(*c) <= 57)
        {
            op_len += *c;
        }

        if (*c == 'S' || *c == 'I')
            op_len = "";

        if (*c == 'M' || *c == 'D')
        {
            referenceSpan += stoi(op_len);
            op_len = "";
        }
    }

    return referenceSpan;
}

/*
 * Initializes variables for read depth calculation
 */
int _readDepth_init( long long unsigned int *totalReadLengths, vector<short> *depthVector, const int currentChr, const bam_hdr_t *bamHdr )
{
    int chromosomeLength;

    *totalReadLengths = 0;
    (*depthVector).clear();

    string refSeqName;
    string targetName = bamHdr->target_name[0];
    if (targetName.find("chr") != string::npos)
    {
        refSeqName = "chr" + to_string(currentChr);
    }
    else
    {
	refSeqName = to_string(currentChr);
    }

    for (int hi = 0; hi < bamHdr->n_targets; hi++)
    {
        if (bamHdr->target_name[hi] == refSeqName)
        {
            chromosomeLength = bamHdr->target_len[hi];
            break;
        }
    }

    for (int i = 0; i < chromosomeLength; i++) {
        (*depthVector).push_back(0);
    }

    return chromosomeLength;
}

/*
 * Adds the number of all matching bases of the read with the reference to the depth vector.
 * Also increments the global read lengths by number of matches. Returns number of matches
 * on success. If the read is in satellite region, it is not included in depth and
 * addRdToDepth returns -1.
 */
int addRdToDepth( const Read *rd, vector<short> *depthVector, unsigned long long int *totalReadLengths )
{
    int op, op_len;
    int matchCount = 0;
    int j = (*rd).pos;

    for (int i = 0; i < (*rd).ncigar; i++)
    {
        op = bam_cigar_op((*rd).cigar[i]);
        if ((op != BAM_CMATCH) && (op != BAM_CDEL)){ continue; }

        op_len = bam_cigar_oplen((*rd).cigar[i]);
        if (op == BAM_CMATCH)
        {
            matchCount += op_len;
            for (int k = 0; k < op_len; k++)
            {
                (*depthVector)[j+k]++;
            }
        }

        j += op_len;
    }

    *totalReadLengths += matchCount;

    return matchCount;
}

/*
 * Finds all alternative mappings of a read. Pushes alternative mappings in the current
 * chromosome to a vector, writes alternative mappings in other chromosomes to a file.
 * Returns 1 on success, 0 otherwise.
 */
int findAltMappings (const int currentChr, string readID, string xa, vector <vector<string> > *altMappings, ofstream *alts )
{
    vector <string> alternatives, fields;
    int chrMismatch;

    splitStringOmitLast( xa, ';', &alternatives );
    alternatives[0].erase(0,1);

    for( int i = 0; i < alternatives.size(); i++ )
    {
        splitStringOmitLast( alternatives[i], ',', &fields );

        chrMismatch = 1;
        try {
            if( stoi(fields[0]) == currentChr )
            {
                chrMismatch = 0;
                (*altMappings).push_back(fields);
            }
        }
        catch( std::exception const & e )
        {
            if( (fields[0] == "MT" &&  currentChr == 25) || (fields[0] == "X" && currentChr == 23) || (fields[0] == "Y" && currentChr == 24) )
            {
                chrMismatch = 0;
                (*altMappings).push_back(fields);
            }
        }

        if( chrMismatch )
        {
            // Write alternative mappings in different chromosomes to an external file
            if ((*alts).is_open())
            {
                (*alts) << readID << "\t";
                for (int i = 0; i < fields.size(); i++)
                {
                    (*alts) << fields[i] << "\t";
                }
                (*alts) << "\n";
            }
            else
                return 0;
        }

        fields.clear();
    }

    return 1;
}
