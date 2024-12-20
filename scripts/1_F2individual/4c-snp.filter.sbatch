#!/bin/bash -l
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-23:00:00
#SBATCH --cpus-per-task=24
module load Java/1.8.0_192
out=2o-MT2B73Ref4
# generate a bam list
vcf=./RX_*/RX.vcf
ref=$refdatabase/GATK/B73Ref4.fa

$path_to_gatk/gatk-4.1.0.0/gatk SelectVariants \
        -R $ref \
        -V $vcf \
        -select 'QUAL>300 && QD > 2.0 && FS < 60.0 && MQ > 40.0 && MQRankSum > -12.5 && ReadPosRankSum > -8.0 && DP >= 100 && DP <= 5000 && AF >= 0.1 && AF <= 0.9 ' \
	--restrict-alleles-to BIALLELIC \
        -select-type SNP \
        -O $out.select.vcf &> $out.select.log

perl ../sourcecode/vcf2AC.pl $out.select.vcf > $out.select.txt

#-select "AF >= 0.1 && AF <= 0.9 && QUAL >= 30.0 && DP >= 100 && DP <= 5000" \
