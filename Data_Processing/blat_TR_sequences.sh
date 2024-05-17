
for f in Ancistrocerus_nigricornis Andrena_haemorrhoa Andrena_minutula Bombus_sylvestris Macropis_europaea Pemphredon_lugubris Vespula_vulgaris Andrena_dorsata Andrena_hattorfiana Apis_mellifera Bombus_terrestris Nysson_spinosus Vespa_crabro
do
	bedtools getfasta -s -fi Genome/$f.fna -bed TR_bed/$f.bed > TR_sequences/$f.fa
done

