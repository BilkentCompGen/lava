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

	-p, --pacbio_bam	: Path to input PacBio BAM file, sorted by coordinate and indexed.
	-i, --ilmn_bam		: Path to input Illumina BAM file, sorted by coordinate and indexed.
	-r, --reference		: 

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

	./lava -p <pacbio.bam> -i <ilmn.bam> -r <GRCh38.fasta> -o <some/path/> -l 1000

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
