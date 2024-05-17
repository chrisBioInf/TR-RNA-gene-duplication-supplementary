while IFS= read -r line; do printf %sn "$line" && fasterq-dump --progress --split-files --threads 8 $line ; done < SraAccList.csv
