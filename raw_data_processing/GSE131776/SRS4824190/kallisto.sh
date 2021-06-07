mkfifo R1.gz R2.gz
curl -Ls ftp://{ftp.sra.ebi.ac.uk/vol1/fastq/SRR913/003/SRR9130163/SRR9130163_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR913/004/SRR9130164/SRR9130164_1.fastq.gz} > R1.gz &
curl -Ls ftp://{ftp.sra.ebi.ac.uk/vol1/fastq/SRR913/003/SRR9130163/SRR9130163_2.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/SRR913/004/SRR9130164/SRR9130164_2.fastq.gz} > R2.gz &
kallisto bus -i /home/Mus_musculus.GRCm38.cdna.all.idx -x 10xv2 -t 4 -o bus_out/ R1.gz R2.gz 
mkdir bus_out/tmp
bustools correct -w /files/10xv2_whitelist.txt -o bus_out/correct_output.bus bus_out/output.bus
bustools sort -t 4 -T bus_out/tmp/ -p bus_out/correct_output.bus | bustools count -o bus_out/genes -g /home/transcripts_to_genes_v2.txt -e bus_out/matrix.ec -t bus_out/transcripts.txt --genecounts -