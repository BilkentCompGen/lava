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

	-p, --pacbio-bam	: Input PacBio BAM file, sorted by coordinate and indexed.
	-i, --ilmn-bam		: Input Illumina BAM file, sorted by coordinate and indexed.
	-r, --reference		: Name of the reference in any format. (e.g. hg19)

Optional:

	-o, --output		: Path to the output bed files. (Default: current dir)
	-l, --min-del-len	: Minimum length for deletions. (Default: 5000)
	-n, --min-inv-len	: Minimum length for inversions. (Default: 10000)

Advanced:

	-c, --min-del-cl-size	: Minimum cluster size to consider for deletions. (Default: 3)
	-z, --min-inv-cl-size	: Minimum cluster size to consider for inversions. (Default: 20)
	-s, --min-del-support	: Minimum number of Illumina support for deletions. (Default: 15)
	-t, --min-inv-support	: Minimum number of Illumina support for inversions. (Default: 70)
	-q, --pacbio-min-qual	: Minimum quality score to consider a PacBio read. (Default: 0)
	-u, --ilmn-min-qual	: Minimum quality score to consider an Illumina read. (Default: 0)

For more information:

	-v, --version		: Display LaVa version.

Example:

	./lava -p pacbio.bam -i ilmn.bam -r b38 -o some/path/ -l 1000

Prebuilt SONIC files (annotations container)
==================================

SONIC files are available under https://github.com/BilkentCompGen/sonic-prebuilt/

 * human_g1k_v37.sonic: SONIC file for Human Reference Genome GRCh37 (1000 Genomes Project version)

 * ucsc_hg19.sonic: SONIC file for the human reference genome, UCSC version build hg19.

 * GRCh38.sonic: SONIC file for the human reference genome build 38.
	
Make sure that the same reference was used to align the reads beforehand (BAM file) and to create the SONIC file. The SONIC files and the reference FASTA files linked above are compatible.

Building a new SONIC file
=======================

Please refer to the SONIC development repository: https://github.com/calkan/sonic/

OUTPUT FORMAT
=============

	chromosome  sv_start_bp1  sv_start_bp2  chromosome  sv_end_bp1  sv_end_bp2  sv_type
