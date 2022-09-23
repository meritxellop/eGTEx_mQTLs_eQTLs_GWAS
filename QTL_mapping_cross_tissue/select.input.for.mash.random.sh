#PBS -N perform_coloc_fastqtl
#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=1GB
#PBS -d /gpfs/data/gtex-group/mpavia/methylation/mashr
#PBS -o /scratch/mpavia/fastQTL/logs/mash_$id_\$PBS_JOBID.out
#PBS -e /scratch/mpavia/fastQTL/logs/mash_$id_\$PBS_JOBID.err

wdir=/gpfs/data/gtex-group/mpavia/methylation
zcat /scratch/mpavia/fastQTL/results/Lung/mqtls/NonPermutedResults/tmpFiles/Lung.$id.regular.txt.gz | shuf -n 200 | awk '{print $1":"$2}' | sort > $wdir/mashr/shuf.$id.txt
join -a 1 $wdir/mashr/shuf.$id.txt  <(zcat /scratch/mpavia/fastQTL/results/BreastMammaryTissue/mqtls/NonPermutedResults/tmpFiles/BreastMammaryTissue.$id.regular.txt.gz | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | sort -k1,1) > $wdir/mashr/BreastMammaryTissue.mQTL_top.fdr005.$id.random.txt
join -a 1 $wdir/mashr/shuf.$id.txt  <(zcat /scratch/mpavia/fastQTL/results/ColonTransverse/mqtls/NonPermutedResults/tmpFiles/ColonTransverse.$id.regular.txt.gz | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | sort -k1,1) > $wdir/mashr/ColonTransverse.mQTL_top.fdr005.$id.random.txt
join -a 1 $wdir/mashr/shuf.$id.txt  <(zcat /scratch/mpavia/fastQTL/results/KidneyCortex/mqtls/NonPermutedResults/tmpFiles/KidneyCortex.$id.regular.txt.gz | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | sort -k1,1) > $wdir/mashr/KidneyCortex.mQTL_top.fdr005.$id.random.txt
join -a 1 $wdir/mashr/shuf.$id.txt  <(zcat /scratch/mpavia/fastQTL/results/Lung/mqtls/NonPermutedResults/tmpFiles/Lung.$id.regular.txt.gz | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | sort -k1,1) > $wdir/mashr/Lung.mQTL_top.fdr005.$id.random.txt
join -a 1 $wdir/mashr/shuf.$id.txt  <(zcat /scratch/mpavia/fastQTL/results/MuscleSkeletal/mqtls/NonPermutedResults/tmpFiles/MuscleSkeletal.$id.regular.txt.gz | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | sort -k1,1) > $wdir/mashr/MuscleSkeletal.mQTL_top.fdr005.$id.random.txt
join -a 1 $wdir/mashr/shuf.$id.txt  <(zcat /scratch/mpavia/fastQTL/results/Testis/mqtls/NonPermutedResults/tmpFiles/Testis.$id.regular.txt.gz | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | sort -k1,1) > $wdir/mashr/Testis.mQTL_top.fdr005.$id.random.txt
join -a 1 $wdir/mashr/shuf.$id.txt  <(zcat /scratch/mpavia/fastQTL/results/Ovary/mqtls/NonPermutedResults/tmpFiles/Ovary.$id.regular.txt.gz | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | sort -k1,1) > $wdir/mashr/Ovary.mQTL_top.fdr005.$id.random.txt
join -a 1 $wdir/mashr/shuf.$id.txt  <(zcat /scratch/mpavia/fastQTL/results/Prostate/mqtls/NonPermutedResults/tmpFiles/Prostate.$id.regular.txt.gz | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | sort -k1,1) > $wdir/mashr/Prostate.mQTL_top.fdr005.$id.random.txt
join -a 1 $wdir/mashr/shuf.$id.txt  <(zcat /scratch/mpavia/fastQTL/results/WholeBlood/mqtls/NonPermutedResults/tmpFiles/WholeBlood.$id.regular.txt.gz | awk '{print $1":"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | sort -k1,1) > $wdir/mashr/WholeBlood.mQTL_top.fdr005.$id.random.txt
