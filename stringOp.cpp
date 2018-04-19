//
// Created by Ezgi Ebren on 22.03.2018.
//

#include "stringOp.h"

/*
 * strtok for c++ strings
 */
void splitStringIncludeLast( string string_to_split, char delimiter, vector <string> *substrings )
{
    string substr;
    size_t pos1, pos2, len;

    pos1 = -1;
    do{
        pos2 = string_to_split.find_first_of( delimiter, pos1+1 );
        len = pos2 - pos1 - 1;
        substr = string_to_split.substr( pos1+1, len );
        (*substrings).push_back(substr);
        pos1 = pos2;

    } while( pos2 != string::npos );
}

/*
 * strtok for c++ strings
 */
void splitStringOmitLast( string string_to_split, char delimiter, vector <string> *substrings )
{
    string substr;
    size_t pos1, pos2, len;

    pos1 = -1;
    pos2 = string_to_split.find_first_of( delimiter, pos1+1 );
    while( pos2 != string::npos )
    {
        len = pos2 - pos1 - 1;
//	cout << "len: " << len << ", ";

        substr = string_to_split.substr( pos1+1, len );
//	cout << "substr: " << substr << endl;
        (*substrings).push_back(substr);

        pos1 = pos2;
        pos2 = string_to_split.find_first_of( delimiter, pos1+1 );
//	cout << "pos2: " << pos2 << ", ";
    }
}

/**
 ** Stores bed file contents in a bed_entry vector.
 ** Returns vector<struct bed_entry>.
 **/
vector<struct bed_entry> file2BedEntry( char *bedFile )
{
    /*
     * Field 0: chr
     * Field 1: start
     * Field 2: end
     * Field 3: type
     */

    int temp, cur;
    string chrPrev;

    vector<string> fields;

    bed_entry entry;
    vector<struct bed_entry> bedContents;

    ifstream f_bed (bedFile);

    cur = 24;
    chrPrev = "-1";
    string line;
    if( f_bed.is_open() )
    {
        while( getline(f_bed, line) )
        {
            splitStringIncludeLast(line, '\t', &fields);

	    try {
		if( fields[0] != chrPrev ) {
		   cur = stoi(fields[0]);
		}
	    }
	    catch( std::exception const & e )
	    {
		if( fields[0] != chrPrev ) {
		    cur++;
		}
	    }
	    entry.chr = cur;

	    try {
		entry.start = stoi(fields[1]);
	    }
	    catch( std::exception const & e )
	    {
		return {};
	    }

	    try {
		entry.end = stoi(fields[2]);
	    }
	    catch( std::exception const & e )
	    {
		return {};
	    }

	    if(  entry.start > entry.end )
	    {
		temp = entry.start;
		entry.start = entry.end;
		entry.end = temp;
	    }

	    if (fields.size() > 3) {
		entry.type = fields[3];
	    }

	    bedContents.push_back(entry);

	    chrPrev = fields[0];
            fields.clear();
        }
    }

    f_bed.close();

    return bedContents;
}
