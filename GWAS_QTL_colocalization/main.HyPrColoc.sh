# Classify GWAS hits into e/mQTL presence categories
# Consider only colocalizations derived from unconditional mQTLs. Hence, from 2734 GWAS hits, consider only 2689 
#R --vanilla < scripts/extract.mQTL.specific.colocs.R

# Select HyPrColoc priors
p1=1e-04
p2=0.98

# For each of the 749 e/mQTL shared hits:
for hit in $(sed 's/\t/,/g' data/mQTL.shared.colocs.GWAS.hits.txt | tail -n+2); do
        gwas=$(echo $hit | awk -F ',' '{print $1}');
        gwas_region=$(echo $hit | awk -F ',' '{print $2}');
# a) Select top-colocalizing (largest PP4) mCpG-tissue case
        top_tissue=$(echo $hit | awk -F ',' '{print $4}');
        top_cpg=$(echo $hit | awk -F ',' '{print $5}');
#        echo mQTL_shared $gwas $gwas_region $top_tissue $top_cpg;
# b) Merge mQTL-GWAS results
# c) For each of the eQTL datasets, merge eQTL-mQTL-GWAS results
# d) Run HyPrColoc on the eQTL(s)-mQTL-GWAS triplet, default parameters, store results 
# e) Repeat operation at a more stringent prior.2
        bash scripts/HyPrColoc.mQTL_specific.cmd $top_tissue $gwas $gwas_region $top_cpg $p1 $p2 mQTL_shared
done


# Hence, instead of 1505, we consider 1460 mQTL specific hits
# For each of the 1460 mQTL specific hits:
for hit in $(sed 's/\t/,/g' data/mQTL.specific.colocs.GWAS.hits.txt | tail -n+2); do
	gwas=$(echo $hit | awk -F ',' '{print $1}');
	gwas_region=$(echo $hit | awk -F ',' '{print $2}');
# a) Select top-colocalizing (largest PP4) mCpG-tissue case
	top_tissue=$(echo $hit | awk -F ',' '{print $4}');
	top_cpg=$(echo $hit | awk -F ',' '{print $5}');
#	echo mQTL_shared $gwas $gwas_region $top_tissue $top_cpg;
# b) Merge mQTL-GWAS results
# c) For each of the eQTL datasets, merge eQTL-mQTL-GWAS results
# d) Run HyPrColoc on the eQTL(s)-mQTL-GWAS triplet, default parameters, store results 
# e) Repeat operation at a more stringent prior.2
	bash scripts/HyPrColoc.mQTL_specific.cmd $top_tissue $gwas $gwas_region $top_cpg $p1 $p2 mQTL_specific
done

# When all jobs completed, parse results.
# When using the variant specific priors, i.e. P∗R=P∗A=0.5, the algorithm will deduce that a set of traits are colocalized when PRxPA≥0.25.

## DISREGARD THIS BLOCK
# Identify as colocalizations, triplets of GWAS-mQTL-eQTL(s) with regional probability (PR) > 0.8 and posterior probability (PP, PP=PRxPA) > 0.5, hence considering PA≥0.625.
for thres in $(echo 0.95 0.98 0.99 0.999); do
	echo $thres; sed "s/\$thres/$thres/g" scripts/extract.external_eQTL.colocs.sh > logs/extract.external_eQTL.colocs.$thres.sh;
	bash logs/extract.external_eQTL.colocs.$thres.sh
	# Identify # GWAS hits corresponding to colocalized triplets of GWAS-mQTL-eQTL(s) clusters
	grep '^mQTL_specific' results/mQTL_specific.summary.$thres.txt | cut -f 1,2,3 | sort | uniq | wc -l
	grep '^mQTL_shared' results/mQTL_shared.summary.$thres.txt | cut -f 1,2,3 | sort | uniq | wc -l
	sed "s/\$thres/$thres/g" scripts/parse.HyPrColoc.results.sh > logs/parse.HyPrColoc.results.$thres.sh
	# Count how many instances of each eQTL dataset exist for mQTL-shared and mQTL-specific colocalizations
	bash logs/parse.HyPrColoc.results.$thres.sh
done

# Run fisher test to check the unequality between mQTL-shared and mQTL-specific colocalizations for over-representation (fisher, on-side, alternative='greater') eQTL sources
# Use prior.2=0.98
R --vanilla < scripts/run.HyPrColoc.Fisher.R

# Classify GWAS hits into e/mQTL presence categories
# Consider only colocalizations derived from unconditional mQTLs. Hence, from 2734 GWAS hits, consider only 2689 
#R --vanilla < scripts/extract.mQTL.specific.colocs.R
## END of DISREGARD THIS BLOCK

# Select HyPrColoc priors
p1=1e-04
p2=0.98

# For each of the 749 e/mQTL shared hits:
for hit in $(sed 's/\t/,/g' data/mQTL.shared.colocs.GWAS.hits.txt | tail -n+2); do
        gwas=$(echo $hit | awk -F ',' '{print $1}');
        gwas_region=$(echo $hit | awk -F ',' '{print $2}');
# a) Select top-colocalizing (largest PP4) mCpG-tissue case
        top_tissue=$(echo $hit | awk -F ',' '{print $4}');
        top_cpg=$(echo $hit | awk -F ',' '{print $5}');
#        echo mQTL_shared $gwas $gwas_region $top_tissue $top_cpg;
# b) Merge mQTL-GWAS results
# c) For each of the eQTL datasets, merge eQTL-mQTL-GWAS results
# d) Run HyPrColoc on the eQTL(s)-mQTL-GWAS triplet, default parameters, store results 
# e) Repeat operation at a more stringent prior.2
	bash scripts/HyPrColoc.mQTL_specific.random_eQTL.cmd $top_tissue $gwas $gwas_region $top_cpg $p1 $p2 mQTL_shared
done


# Hence, instead of 1505, we consider 1460 mQTL specific hits
# For each of the 1460 mQTL specific hits:
for hit in $(sed 's/\t/,/g' data/mQTL.specific.colocs.GWAS.hits.txt | tail -n+2); do
	gwas=$(echo $hit | awk -F ',' '{print $1}');
	gwas_region=$(echo $hit | awk -F ',' '{print $2}');
# a) Select top-colocalizing (largest PP4) mCpG-tissue case
	top_tissue=$(echo $hit | awk -F ',' '{print $4}');
	top_cpg=$(echo $hit | awk -F ',' '{print $5}');
#	echo mQTL_shared $gwas $gwas_region $top_tissue $top_cpg;
# b) Merge mQTL-GWAS results
# c) For each of the eQTL datasets, merge eQTL-mQTL-GWAS results
# d) Run HyPrColoc on the eQTL(s)-mQTL-GWAS triplet, default parameters, store results 
# e) Repeat operation at a more stringent prior.2
	bash scripts/HyPrColoc.mQTL_specific.random_eQTL.cmd $top_tissue $gwas $gwas_region $top_cpg $p1 $p2 mQTL_specific
done

# When all jobs completed, parse results to calculate FPR.
# Annotate FPR
bash run.HyPrColoc.annotate.fpr.sh

# Count hits
bash count.HyPrColoc.hits.sh
