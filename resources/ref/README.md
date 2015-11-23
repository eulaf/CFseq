# CFTR reference sequence

This directory contains the reference sequence for the CFTR gene.

## Files

- `CFTR_genomic_refseq.fasta` - FASTA sequence with newlines removed.
- `CFTR_genomic_refseq.fasta.fai` - FASTA index file required by GATK.  
  Created using samtools' faidx command.
- `CFTR_genomic_refseq.dict` - FASTA sequence dictionary file used by GATK.  
  Created by running Picard's CreateSequenceDictionary function.
- `CFTR_genomic_refseq.txt` - original sequence file
- `CFTR_genomic_refseq.fasta.[sa,pac,bwt,ann,amb]` - bwa index files

## Commands

To create samtools fasta index: 
``` 
samtools faidx CFTR_genomic_refseq.fasta 
```

To create reference dict: 
```
java -jar picard.jar CreateSequenceDictionary R=CFTR_genomic_refseq.fasta O=CFTR_genomic_refseq.dict
```

To create bwa index files: 
```
bwa index CFTR_genomic_refseq.fasta
```

