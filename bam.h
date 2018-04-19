#ifndef UNTITLED_BAM_H
#define UNTITLED_BAM_H

#include <getopt.h>

#include "stringOp.h"
#include "read.h"
#include "sv.h"
#include "cluster.h"

using namespace std;

int readIlluminaBam( const samFile *fp_in, bam_hdr_t *bamHdrIl, bam1_t *aln, const hts_idx_t *bamIndex, sonic *snc, const int currentChr, const int chromosomeLength,
                      vector<SV> *deletions, vector<SV> *inversions );

int readBam( vector<SV> *deletions_final, vector<SV> *inversions_final, char *sonic_filename );

#endif //UNTITLED_BAM_H
