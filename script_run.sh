nohup cat index.txt | while read id; do bash bash_script/split_script/split_script_${id} &>logs/${id}.log & done
# for multithreading accelerating, 83 threads