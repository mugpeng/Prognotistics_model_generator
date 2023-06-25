seq 1 1763 | while read id; do echo "Rscript main2.R $id"; done &> all_script.sh
# 9034 is the length of selected msigdb gene pathways.