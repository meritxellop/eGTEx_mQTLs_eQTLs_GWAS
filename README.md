[![DOI](https://zenodo.org/badge/539703495.svg)](https://zenodo.org/badge/latestdoi/539703495)
# 1. Introduction
Tissue-specific characterization of DNA methylation (DNAm) is needed to understand its role in gene regulation and its relevance to complex traits. We generated array-based DNAm profiles for 987 human samples from the Genotype-Tissue Expression (GTEx) project, representing 9 tissue types and 424 subjects. We characterized methylome and transcriptome correlations (eQTMs), genetic regulation in cis (mQTLs and eQTLs) across tissues, and e/mQTLs links to complex traits. With this DNAm-focused integrative analysis, we contribute to the understanding of molecular regulatory mechanisms in human tissues and their impact on complex traits.


# 2. Contents

## DNAm_processing
Includes code for DNAm processing with ChAMP(v.2.8.6) and minfi(1.36.0), available at Bioconductor https://www.bioconductor.org/.  

## QTL_mapping_single_tissue
Includes code for primary QTL detection using FastQTL(v.2.184). The adaptation of FastQTL is available at https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl. Scripts for conditional QTL mapping have been adapted from https://github.com/funpopgen/multiple_eqtl_mapping scripts.

## QTL_mapping_cross_tissue
Includes code for cross-tissue QTL mapping by leveraging QTL signal across tissues with mashr. mashr is available at https://github.com/stephenslab/mashr.  

## TORUS_enrichment
Includes code for QTL-annotation enrichment analysis using torus. torus is available at https://github.com/xqwen/torus.

## eQTM_mapping_single_and_cross_tissue
Includes code for single-tissue and cross-tissue eQTM mapping by leveraging eQTM signal across tissues with mashr. mashr is available at https://github.com/stephenslab/mashr.

## GWAS_QTL_colocalization
Includes code for QTL-GWAS colocalization analysis using coloc, fastenloc and hyprcoloc. fastenloc is available at https://github.com/xqwen/fastenloc. coloc is available at https://github.com/chr1swallace/coloc/. hyprcoloc is available at https://github.com/jrs95/hyprcoloc.

## Figures
Includes code to generate manuscript figures. Data can be found at https://figshare.com/projects/DNA_methylation_QTL_mapping_across_diverse_human_tissues_provides_molecular_links_between_genetic_variation_and_complex_traits/149524

