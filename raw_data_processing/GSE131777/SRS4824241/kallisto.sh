rm *.gz
mkfifo barcode.gz bio.gz
fastq-dump --split-files --gzip SRR9130263
fastq-dump --split-files --gzip SRR9130264
fastq-dump --split-files --gzip SRR9130265
fastq-dump --split-files --gzip SRR9130266
cat *1.fastq.gz > barcode.gz &
cat *2.fastq.gz > bio.gz &
kallisto bus -i /scratch/mfiruleva/winter/genomes/mus/Mus_musculus.GRCm38.cdna.all.idx -x 10xv2 -t 4 -o bus_out/ barcode.gz bio.gz
mkdir bus_out/tmp
bustools correct -w /scratch/mfiruleva/GSEs_processing/wget_and_process/10xv2_whitelist.txt -o bus_out/correct_output.bus bus_out/output.bus
rm bus_out/output.bus
bustools sort -t 4 -T bus_out/tmp/ -p bus_out/correct_output.bus | bustools count -o bus_out/genes -g /scratch/mfiruleva/winter/genomes/mus/transcripts_to_genes_v2.txt -e bus_out/matrix.ec -t bus_out/transcripts.txt --genecounts -