for tissue in $(grep -v Kidney ../coloc/data/tissues.txt); do 
	TISSUE=$(echo $tissue | sed 's/_//g'); 
	echo $tissue; 
	sed "s/TISSUE/$TISSUE/g" scripts/summary.eQTM.results.cmd  > logs/summary.eQTM.results.$tissue.cmd;
	if [[ ! -s results/$TISSUE.cor.stats.complete.txt ]]; then 
		echo TODO $tissue; 
		qsub logs/summary.eQTM.results.$tissue.cmd;
		ls logs/summary.eQTM.results.$tissue.cmd;
	fi;
done;
