perl ./bin/Data_Prepare.pl -I ./data/common/hg19.fa -P CCGG -S 1 -O ./data/common/hg19.fragment.bed -B ./bin
./bin/bismark_genome_preparation  --bowtie2   --path_to_bowtie ./bin/  --verbose ./data/common/
