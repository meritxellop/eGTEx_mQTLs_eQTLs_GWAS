cat <(head -n 1 output/TF/meta.summary.TF.full.txt) <(grep -f <(cat <(sort -gk11  output/TF/meta.summary.TF.full.txt | grep eQTL | awk '{if($3>0 && $11<0.01) {print}}' | sort -rgk3 | head -n 15 | cut -f 2) <(sort -gk11 output/TF/meta.summary.TF.full.txt | grep mQTL | awk '{if($3>0 && $11<0.01) {print}}' | sort -rgk3 | head -n 15 | cut -f 2) | sort | uniq) output/TF/meta.summary.TF.full.txt)
