seq 1 11869 | while read id; do echo "Rscript main.R $id"; done &> all_script.sh
# 11869 is the length of selected msigdb gene pathways.