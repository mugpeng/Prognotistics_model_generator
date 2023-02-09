nohup cat index.txt | while read id; do bash bash_script/split_script/split_script_${id} &>logs/${id}.log & done
# for multithreading accelerating, 83 threads

# nohup cat index.txt | head -36 | while read id; do bash bash_script/split_script/split_script_${id} &>logs/${id}.log & done
nohup cat index.txt | sed '1,66d' | while read id; do bash bash_script/split_script/split_script_${id} &>logs/${id}.log & done