# CFseq
CFseq analysis pipeline scripts.

Lefterova MI, Shen P, Odegaard JI, et al.
Next-generation molecular testing of newborn dried blood spots for cystic fibrosis.
*J Mol Diagn* (http://jmd.amjpathol.org).


## Required software

The pipeline uses the following software and was designed using the version numbers listed.  
All software are freely available for non-commercial use.

- ANNOVAR:
    * http://www.openbioinformatics.org/annovar/
	* release 150322

- BWA:
    * http://bio-bwa.sourceforge.net
	* version 0.7.12

- FreeBayes: 
    * https://github.com/ekg/freebayes
	* version v0.9.21-7-g7dd41db

- GATK: 
    * http://www.broadinstitute.org/gatk/ 
	* version 3.3

- Picard
    * http://broadinstitute.github.io/picard/
	* version 1.124

- SAMtools:
    * http://www.htslib.org
	* version 1.1

- cross_match/phrap:
    * http://www.phrap.org/phredphrapconsed.html
	* release 990329

## Setup

### bin:

Place these executables in the bin directory:

* annotate_variation.pl - from ANNOVAR
* convert2annovar.pl - from ANNOVAR
* table_annovar.pl - from ANNOVAR
* bwa
* cross_match.manyreads
* freebayes
* GenomeAnalysisTK.jar
* picard.jar
* samtools

### resources/annovar_dbs:

Download the 'refGene' database from ANNOVAR to resources/annovar_dbs.


