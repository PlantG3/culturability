bamdir=../3-aln.filter2B73Ref4
ref=$refdatabase/B73Ref4/GATK/B73Ref4.fa

perl ../soucecode/gatk.sbatch.pl --outbase RX \
  --bampaths $bamdir \
  --ref $ref \
  --mem 12G --maxlen 1000000
  #--checkscript
