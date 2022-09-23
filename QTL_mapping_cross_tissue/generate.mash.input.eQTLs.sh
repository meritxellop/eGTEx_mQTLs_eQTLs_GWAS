file=pvalue.eQTL_top.fdr005.all.txt; echo -e "BreastMammaryTissue\tColonTransverse\tKidneyCortex\tLung\tMuscleSkeletal\tOvary\tProstate\tTestis\tWholeBlood" > $file;paste <(awk '{print $1,$6}' BreastMammaryTissue.eQTL_top.fdr005.all.txt | tr ' ' '\t' | sed "s/\t$/\tNA/g") <(awk '{print $6}' ColonTransverse.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $6}' KidneyCortex.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $6}' Lung.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $6}' MuscleSkeletal.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $6}' Ovary.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $6}' Prostate.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $6}' Testis.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $6}' WholeBlood.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") >> pvalue.eQTL_top.fdr005.all.txt

file=beta.eQTL_top.fdr005.all.txt; echo -e "BreastMammaryTissue\tColonTransverse\tKidneyCortex\tLung\tMuscleSkeletal\tOvary\tProstate\tTestis\tWholeBlood" > $file;paste <(awk '{print $1,$7}' BreastMammaryTissue.eQTL_top.fdr005.all.txt | tr ' ' '\t' | sed "s/\t$/\tNA/g") <(awk '{print $7}' ColonTransverse.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $7}' KidneyCortex.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $7}' Lung.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $7}' MuscleSkeletal.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $7}' Ovary.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $7}' Prostate.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $7}' Testis.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $7}' WholeBlood.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") >> beta.eQTL_top.fdr005.all.txt

file=beta_se.eQTL_top.fdr005.all.txt; echo -e "BreastMammaryTissue\tColonTransverse\tKidneyCortex\tLung\tMuscleSkeletal\tOvary\tProstate\tTestis\tWholeBlood" > $file;paste <(awk '{print $1,$8}' BreastMammaryTissue.eQTL_top.fdr005.all.txt | tr ' ' '\t' | sed "s/\t$/\tNA/g") <(awk '{print $8}' ColonTransverse.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $8}' KidneyCortex.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $8}' Lung.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $8}' MuscleSkeletal.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $8}' Ovary.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $8}' Prostate.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $8}' Testis.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") <(awk '{print $8}' WholeBlood.eQTL_top.fdr005.all.txt | sed "s/^$/NA/g") >> beta_se.eQTL_top.fdr005.all.txt
