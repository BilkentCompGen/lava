/*
************************************************************
**
** Author: Ezgi Ebren - 2018
**
**
** This program detects and reports large deletions and
** inversions.
**
************************************************************
**
** Usage:
**
** ./lava -p <pacbio_bam> -i <illumina_bam> -r <name_of_reference_sequence>
**
**/


/*
 * Finds deletions in Illumina data using read-pairs.
 */

#include "bam.h"

int readIlluminaBam( const samFile *fp_in, bam_hdr_t *bamHdrIl, bam1_t *aln, const hts_idx_t *bamIndex, sonic *snc, const int currentChr, const int chromosomeLength,
                      vector<SV> *deletions, vector<SV> *inversions )
{
    vector<SV> forwardPairs, reversePairs;

    Params *params;

    int unmapped, mateUnmapped, inReverseStrand, mateInReverseStrand, notPrimary, qcFail, dup, suppl;

    params = Params::getInstance();

    hts_itr_t *iter = bam_itr_queryi(bamIndex, currentChr-1, 0, chromosomeLength);

    while( bam_itr_next(fp_in, iter, aln) > 0)
    {
        if (aln->core.tid+1 != currentChr)
        {
            return -1;
        }

        unmapped = aln->core.flag&BAM_FUNMAP;
        mateUnmapped = aln->core.flag&BAM_FMUNMAP;
        inReverseStrand = aln->core.flag&BAM_FREVERSE;
        mateInReverseStrand = aln->core.flag&BAM_FMREVERSE;
        notPrimary = aln->core.flag&BAM_FSECONDARY;
        qcFail = aln->core.flag&BAM_FQCFAIL;
        dup = aln->core.flag&BAM_FDUP;
        suppl = aln->core.flag&BAM_FSUPPLEMENTARY;
        if ( unmapped || mateUnmapped ) //|| notPrimary || qcFail || dup || suppl )
        {
            continue;
        }

        if ( aln->core.qual >= params->ilmn_min_qual )
        {
            if ( (inReverseStrand == 0) && (mateInReverseStrand != 0) && (aln->core.mpos - aln->core.pos > 0) )
            {
                findIlluminaDel(deletions, aln, currentChr, snc);
            }
        }

        if ( aln->core.qual >= params->ilmn_min_qual )
        {
            if( ( ( (inReverseStrand == 0) && (mateInReverseStrand == 0) ) || ( (inReverseStrand != 0) && (mateInReverseStrand != 0) ) ) && ( aln->core.mpos - aln->core.pos > 0 ) )
            {
                findIlluminaInv(inversions, &forwardPairs, &reversePairs, aln, currentChr, snc);
            }
        }
    }
}

/**
 */
int readBam( vector<SV> *deletions_final, vector<SV> *inversions_final, char *sonic_filename )
{
    samFile *fp_pb_in, *fp_il_in;
    bam_hdr_t *bamHdrPb, *bamHdrIl;
    bam1_t alignment, *aln_pb, *aln_il;
    hts_idx_t *ilBamIndex;

    Read rd;
    SV sv;

    Unmmap candidates;

    vector<SV> deletions_pb, deletions_il;
    vector<SV> inversions_pb, inversions_il;
    vector<Cluster> clusters_del, clusters_inv;

    forward_list<uint32_t *> cigarPb;

    int unmapped, notPrimary, qcFail, dup, suppl;
    uint64_t start, len, tol;

    vector<short> depthVector;
    int currentChr, chromosomeLength;
    int read_length;
    unsigned long long int totalReadLengths;
    double chrReadDepth;
    string refSeqName;

    vector<int> SClengths;

    string finalDelsBed;
    string finalInvsBed;

    sonic *snc = sonic_load(sonic_filename);

    Params *params = Params::getInstance();

    if (params->output == "") {
	finalDelsBed = "lava_deletions.bed";
	finalInvsBed = "lava_inversions.bed";
    }
    else {
	finalDelsBed = params->output + "_lava_deletions.bed";
	finalInvsBed = params->output + "_lava_inversions.bed";
    }

    ofstream f_dels_all (finalDelsBed);
    ofstream f_invs_all (finalInvsBed);

    fp_pb_in = hts_open((params->pacbio_bam_file).c_str(), "r"); //open pacbio bam file
    bamHdrPb = sam_hdr_read(fp_pb_in); //read header
    aln_pb = bam_init1(); //initialize an alignment

    fp_il_in = hts_open((params->ilmn_bam_file).c_str(), "r"); //open illumina bam file
    bamHdrIl = sam_hdr_read(fp_il_in); //read header
    aln_il = bam_init1(); //initialize an alignment
    ilBamIndex = sam_index_load(fp_il_in, (params->ilmn_bam_file).c_str());

    int satellite = 0;

    currentChr = -1;
    int sam_flag = sam_read1(fp_pb_in, bamHdrPb, aln_pb);
    while (sam_flag > 0)
    {
        if (currentChr == -1)
        {
            currentChr = aln_pb->core.tid+1;
            chromosomeLength = _readDepth_init(&totalReadLengths, &depthVector, currentChr, bamHdrPb);
        }

        if(currentChr == aln_pb->core.tid+1)
        {
            rd.flag = aln_pb->core.flag;

            unmapped = rd.flag & BAM_FUNMAP;
            notPrimary = rd.flag & BAM_FSECONDARY;
            qcFail = rd.flag & BAM_FQCFAIL;
            dup = rd.flag & BAM_FDUP;
            suppl = rd.flag & BAM_FSUPPLEMENTARY;
            if (unmapped != 0 || aln_pb->core.pos <= 0)// || notPrimary || qcFail || dup || suppl )
            {
		sam_flag = sam_read1(fp_pb_in, bamHdrPb, aln_pb);
                continue;
            }
            rd.qname = (std::string) bam_get_qname(aln_pb);  // Read ID
            rd.rname = currentChr;                           // Chromosome
            rd.pos = aln_pb->core.pos;
            rd.ncigar = aln_pb->core.n_cigar;

            uint32_t *data = new uint32_t[rd.ncigar];

            //Copies the cigar section of alignment to an array
            //and stores it in a forward list
	    int i = 0;
	    for (uint32_t *p = (uint32_t *)(aln_pb->data + aln_pb->core.l_qname); i < (rd.ncigar); p++)
            {
                data[i] = *p;
                i++;
            }

            cigarPb.emplace_front(data);
            rd.cigar = *(cigarPb.begin());

            read_length = readLength(&rd);

            if (sonic_is_satellite(snc, snc->chromosome_names[currentChr-1], rd.pos, rd.pos+read_length) == 0) {
		addRdToDepth(&rd, &depthVector, &totalReadLengths);
            }
	    else{
		satellite++;
	    }

            if( params->min_del_len <= 2000 && delCigar(&rd, &start, &len, &tol) == 1 )
            {
                int end = start + len - 1;

                sv.chromosome = rd.rname + 1;
                sv.pos.start1 = start;
                sv.pos.start2 = start;
                sv.pos.end1 = end;
                sv.pos.end2 = end;
                sv.tolerance = tol;
                sv.readID = rd.qname;
                sv.type = deletion;

                deletions_pb.push_back(sv);

                // continue;
            }

            // Discard reads that are < 300 bp in length, since they may be ALU
            if (read_length < 300) {
		sam_flag = sam_read1(fp_pb_in, bamHdrPb, aln_pb);
		continue;
	    }

            softClipLengths(&rd, &SClengths);
            if (SClengths[0] == 0 && SClengths[1] == 0) {
		sam_flag = sam_read1(fp_pb_in, bamHdrPb, aln_pb);
		continue;
	    }

            //Inserts the read in an std::unordered_multimap.
            //After this step, all candidate reads will be present in the multimap.
            candidates.insert(Unmmap::value_type(rd.qname, rd));
        }
        else
        {
            //Calculate avg read depth of the chromosome
            chrReadDepth = ((double) (totalReadLengths)) / (2.0 * (double) (chromosomeLength));
            findPacbioSVs (&candidates, &deletions_pb, &inversions_pb, currentChr, snc);
	    readIlluminaBam (fp_il_in, bamHdrIl, aln_il, ilBamIndex, snc, currentChr, chromosomeLength, &deletions_il, &inversions_il);
	    clusterSVs(&deletions_pb, &clusters_del);
	    clusterSVs(&inversions_pb, &clusters_inv);
	    filterSVs(&deletions_il, &clusters_del, deletions_final, &depthVector, chrReadDepth, currentChr, &f_dels_all);
	    filterSVs(&inversions_il, &clusters_inv, inversions_final, &depthVector, chrReadDepth, currentChr, &f_invs_all);

            for (auto it = cigarPb.begin(); it != cigarPb.end(); ++it)
            {
                delete[] (*it);
            }

	    candidates.clear();
	    deletions_pb.clear();
	    deletions_il.clear();
	    inversions_pb.clear();
	    inversions_il.clear();
	    clusters_del.clear();
	    clusters_inv.clear();

	    cigarPb.clear();

	    depthVector.clear();
	    SClengths.clear();

	    /*** START READING THE NEXT CHROMOSOME ***/
	    currentChr = aln_pb->core.tid + 1;

	    if (currentChr > 24)
	    {
		return 1;
	    }

            chromosomeLength = _readDepth_init(&totalReadLengths, &depthVector, currentChr, bamHdrPb);
            continue;
        }

        sam_flag = sam_read1(fp_pb_in, bamHdrPb, aln_pb);
    }

    //Calculate avg read depth of the chromosome
    chrReadDepth = ((double) (totalReadLengths)) / (2.0 * (double) (chromosomeLength));
    findPacbioSVs (&candidates, &deletions_pb, &inversions_pb, currentChr, snc);
    readIlluminaBam (fp_il_in, bamHdrIl, aln_il, ilBamIndex, snc, currentChr, chromosomeLength, &deletions_il, &inversions_il);
    clusterSVs(&deletions_pb, &clusters_del);
    clusterSVs(&inversions_pb, &clusters_inv);
    filterSVs(&deletions_il, &clusters_del, deletions_final, &depthVector, chrReadDepth, currentChr, &f_dels_all);
    filterSVs(&inversions_il, &clusters_inv, inversions_final, &depthVector, chrReadDepth, currentChr, &f_invs_all);

    //Close and destroy objects

    for (auto it = cigarPb.begin(); it != cigarPb.end(); ++it)
    {
        delete[] (*it);
    }

    sam_close(fp_pb_in);
    bam_destroy1(aln_pb);
    bam_hdr_destroy(bamHdrPb);

    sam_close(fp_il_in);
    bam_destroy1(aln_il);
    bam_hdr_destroy(bamHdrIl);

    candidates.clear();
    deletions_pb.clear();
    deletions_il.clear();
    inversions_pb.clear();
    inversions_il.clear();
    clusters_del.clear();
    clusters_inv.clear();

    cigarPb.clear();

    f_dels_all.close();
    f_invs_all.close();

    cout << "Done." << endl;
}

// Gets command line parameters
int parse_command_line( int argc, char* argv[] )
{
    static struct option long_options[] =
    {
	{"help"			, no_argument,	     0, 'h'},
	{"version"		, no_argument,	     0, 'v'},
	{"pacbio_bam"		, required_argument, 0, 'p'},
	{"ilmn_bam"		, required_argument, 0, 'i'},
	{"reference"		, required_argument, 0, 'r'},
	{"output"		, required_argument, 0, 'o'},
	{"min_del_len"		, required_argument, 0, 'l'},
	{"min_inv_len"		, required_argument, 0, 'n'},
	{"min_del_cluster_size"	, required_argument, 0, 'c'},
	{"min_inv_cluster_size"	, required_argument, 0, 'z'},
	{"min_del_support"	, required_argument, 0, 's'},
	{"min_inv_support"	, required_argument, 0, 't'},
	{"pacbio_min_qual"	, required_argument, 0, 'q'},
	{"ilmn_min_qual"	, required_argument, 0, 'u'},
	{0			, 0		   , 0,	 0 }

    };

    int o, index;
    Params *params;

    if (argc == 1)
    {
//	print_help();
	return -1;
    }

    //TODO check parameter validity

    params = Params::getInstance();
    while ((o = getopt_long(argc, argv, "hvp:i:r:o:l:n:c:z:s:t:q:u:", long_options, &index)) != -1)
    {
	switch(o)
	{
	    case 'p':
		params->pacbio_bam_file = string(optarg);
	    break;
	    case 'i':
		params->ilmn_bam_file = string(optarg);
	    break;
	    case 'r':
		params->name_of_reference = string(optarg);
	    break;
	    case 'o':
		params->output = string(optarg);
	    break;
	    case 'l':
		params->min_del_len = atoi(optarg);
	    break;
	    case 'n':
		params->min_inv_len = atoi(optarg);
	    break;
	    case 'c':
		params->min_del_cluster_size = atoi(optarg);
	    break;
	    case 'z':
		params->min_inv_cluster_size = atoi(optarg);
	    break;
	    case 's':
		params->min_del_support = atoi(optarg);
	    break;
	    case 't':
		params->min_inv_support = atoi(optarg);
	    break;
	    case 'q':
		params->pacbio_min_qual = atoi(optarg);
	    break;
	    case 'u':
		params->ilmn_min_qual = atoi(optarg);
	    break;
	    case 'h':
//		print_help();
		return -1;
	    break;
	    case 'v':
		cout << "LaVa: Large Variation discovery using hybrid sequence data." << endl;
		cout << "Version " << LAVA_VERSION << "\tLatest update: " << LAVA_UPDATE << ", build date: " << BUILD_DATE << endl;
		return -1;
	    break;
	}
    }

    if (params->pacbio_bam_file == "")
    {
	cout << "[LAVA PARAMETER ERROR] Input PacBio BAM file must be given using the --pacbio_bam or -p option." << endl;
	return -1;
    }

    if (params->ilmn_bam_file == "")
    {
	cout << "[LAVA PARAMETER ERROR] Input Illumina BAM file must be given using the --ilmn_bam or -i option." << endl;
    }

    if (params->name_of_reference == "")
    {
	cout << "[LAVA PARAMETER ERROR] Name of the reference sequence must be given using the --reference or -r option." << endl;
	return -1;
    }
}

int main(int argc, char *argv[])
{
    cout << "LaVa: Detection of Large Variations.\n" << endl;
    cout << "Start." << endl;

    vector<SV> deletions_final, inversions_final;

    string sonic_filename;

    Params *params;
    if (!parse_command_line(argc, argv)) {
	return -1;
    }

    params = Params::getInstance();

    if (params->name_of_reference.find("hg") != string::npos && params->name_of_reference.find("19") != string::npos)
    {
	sonic_filename = "aux/ucsc_hg19.sonic";
    }
    else if (params->name_of_reference.find("b37") != string::npos || params->name_of_reference.find("v37") != string::npos)
    {
	sonic_filename = "aux/human_g1k_v37.sonic";
    }
    else if (params->name_of_reference.find("b38") != string::npos || params->name_of_reference.find("hg38") != string::npos)
    {
	 sonic_filename = "aux/ucsc_hg38.sonic";
    }
    else {
	cout << "Please build a new SONIC file for the reference you are using and place it in the directory \"LaVa/aux/ \". More info: https:/" << "/github.com/calkan/sonic/" << endl;
	return -1;
    }

    readBam(&deletions_final, &inversions_final, (char *)sonic_filename.c_str());

    cout << "Done." << endl;
    return 0;
}
