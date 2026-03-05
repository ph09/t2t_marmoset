# Variant burden at Alzheimer's disease loci

Variant effect predictions were done using [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/) version 5.4a (build 2025-11-25 12:22). 
SnpEff build was used to construct a custom genome database for annotation using reference genomic fasta and gene transfer format (GTF v2.2) files for the new T2T assembly, following the instructions described [here](https://pcingola.github.io/SnpEff/snpeff/build_db/). 

```{bash}
# here, caljac240 is the name of the custom genome built by us and added in the snpEff.config file.
java -Xms4g -Xmx48g -jar ./snpEff.jar build -gtf22 -noCheckCds -noCheckProtein -v caljac240 
```

A bed file was constructed by extracting genomic coordinates from 'gene' entries in the .gtf file for [81 genes](./gene_loci_targets.txt) considered based on their association with Alzheimer’s disease and related pathology (Bellenguez et al., 2022, others). 
This bed file was used to extract variants at loci of interest from the joint VCF file using bcftools view (bcftools version 1.9; using htslib 1.9). This subset of variants was further annotated with SnpEff. 

```{bash}
# AD_loci_subset_N230.vcf.gz is the subset of variants filtered for the AD loci of our interest.
java -Xms4g -Xmx48g -jar snpEff.jar -v caljac240 AD_loci_subset_N230.vcf.gz > annotated_AD_loci_N230.vcf.gz
```

SnpEff produces two output files apart from an annotated VCF file: 
1) an HTML report summarizing sequence level changes and annotated variant impacts and effects; and
2) a TSV file showing variant annotations at the gene level.

For the variant burden at AD loci represented in the T2T marmoset paper (Hebbar et al., 2026), the TSV file detailing variant types and impacts per transcript ID was then utilized for analysis, with multiple transcripts collapsed into one entry per gene for representative purposes. 
[This Quarto document](./plotting_variant_burden.qmd) contains the code needed to reproduce the figures in the paper.
