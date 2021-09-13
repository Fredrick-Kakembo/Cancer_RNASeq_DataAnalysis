for sample in `cat sample_list.txt`
do
	fastq-dump --gzip --split-3 $sample
done
