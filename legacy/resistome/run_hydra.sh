samtools view -uF 2 tmp.bam  | bamToFastq -i - -fq sample.tier1.disc.1.fq -fq2 sample.tier1.disc.2.fq
#novoindex  -k 14 -s 1 -m novoindexfile test.fasta
novoalign -c 1 -d novoindexfile -f sample.tier1.disc.1.fq sample.tier1.disc.2.fq -i 500 50 -r Random -o SAM | samtools view -Sb - > sample.tier2.queryorder.bam
samtools view -uF 2 sample.tier2.queryorder.bam | bamToFastq -i - -fq sample.tier2.disc.1.fq -fq2 sample.tier2.disc.2.fq
novoalign -c 1 -d novoindexfile -f sample.tier2.disc.1.fq sample.tier2.disc.2.fq -i 500 50 -r Ex 1100 -t 300 -o SAM | samtools view -Sb - > sample.tier3.queryorder.bam
bamToBed -i sample.tier3.queryorder.bam -tag NM | ~/Hydra-Version-0.5.3/scripts/pairDiscordants.py -i stdin -m hydra -z 800 > sample.disc.bedpe
~/Hydra-Version-0.5.3/scripts/dedupDiscordants.py -i sample.disc.bedpe -s 3 > sample.disc.deduped.bedpe
~/Hydra-Version-0.5.3/bin/hydra -in sample.disc.bedpe -out sample.breaks -mld 500 -mno 1500