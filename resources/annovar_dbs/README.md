## annovar_dbs

Expected files:
- hg19_refGeneMrna.fa
- hg19_refGene.txt
- hg19_refGeneVersion.txt

To download 'refGene' database from ANNOVAR:

```
annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene .
```

