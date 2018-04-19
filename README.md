LAVA
====

Large Variation discovery using hybrid sequence data

Dependencies
=============

 * htslib (http://htslib.org/)
 * sonic  (https://github.com/calkan/sonic)

Compilation
===========

	make libs
	
	make
	

Running LaVa
============

Required:

	-p, --pacbio_bam	: Input PacBio BAM file, sorted by coordinate and indexed.
	-i, --ilmn_bam		: Input Illumina BAM file, sorted by coordinate and indexed.
	-r, --reference		: Name of the reference in any format. (e.g. hg19)

Optional:

	-o, --output		: Path to the output bed files. (Default: current dir)
	-l, --min_del_len	: Minimum length for deletions. (Default: 5000)
	-n, --min_inv_len	: Minimum length for inversions. (Default: 10000)

Advanced:

	-c, --min_del_cl_size	: Minimum cluster size to consider for deletions. (Default: 3)
	-z, --min_inv_cl_size	: Minimum cluster size to consider for inversions. (Default: 20)
	-s, --min_del_support	: Minimum number of Illumina support for deletions. (Default: 15)
	-t, --min_inv_support	: Minimum number of Illumina support for inversions. (Default: 70)
	-q, --pacbio_min_qual	: Minimum quality score to consider a PacBio read. (Default: 0)
	-u, --ilmn_min_qual	: Minimum quality score to consider an Illumina read. (Default: 0)

For more information:

	-v, --version		: Display LaVa version.

Example:

	./lava -p pacbio.bam -i ilmn.bam -r b38 -o some/path/ -l 1000

Building a new SONIC file
=======================

Please refer to the SONIC development repository: https://github.com/calkan/sonic/

OUTPUT FORMAT
=============

chromosome  sv_start_bp1  sv_start_bp2  chromosome  sv_end_bp1  sv_end_bp2  sv_type
