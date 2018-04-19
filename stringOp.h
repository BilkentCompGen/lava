#ifndef UNTITLED_STRINGOP_H
#define UNTITLED_STRINGOP_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <forward_list>
#include <cmath>
#include <iomanip>
#include <typeinfo>
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"

#include "sonic.h"

using namespace std;

struct bed_entry {
    int chr;
    int start;
    int end;
    string type;

    bed_entry () {
	chr = -1;
	start = 0;
	end = 0;
	type = "none";
    }
};

/*
 * strtok for c++ strings
 */
void splitStringIncludeLast( string string_to_split, char delimiter, vector <string> *substrings );

/*
 * strtok for c++ strings
 */
void splitStringOmitLast( string string_to_split, char delimiter, vector <string> *substrings );

/**
 ** Stores bed file contents in a bed_entry vector.
 ** Returns vector<struct bed_entry>.
 **/
vector<struct bed_entry> file2BedEntry( char *bedFile );

#endif //UNTITLED_STRINGOP_H
