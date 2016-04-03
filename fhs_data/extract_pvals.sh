for i in phs*.txt; do awk -F "\t" {print \t } $i > ${i}_pvals; done
for i in *.txt_pvals; do tail -n +23 $i > ${i}_tmp; mv ${i}_tmp $i; done
