#!/bin/bash
star_ref=/homes/liu3zhen/references/B73Ref4/STAR

perl sourcecode/STAR.sbatch.pl \
	--java Java/1.8.0_192 \
	--mem 3 --threads 16 --time 12:00:00 \
	--star_cmd path_to_your/STAR \
	--indir ../2B_trim \
	--dbdir $star_ref \
	--fq1feature .R1.pair.fq \
	--fq2feature .R2.pair.fq \
	--alignIntronMax 50000 \
	--alignMatesGapMax 50000 \
	--outSAMattrIHstart 0 \
	--outSAMmultNmax 1 \
	--outSAMstrandField intronMotif  \
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--outSAMtype "BAM SortedByCoordinate" \
	--quantMode GeneCounts \
	--outFilterMismatchNmax 5 \
	--outFilterMismatchNoverLmax 0.05 \
	--outFilterMatchNmin 80 \
	--outSJfilterReads Unique \
	--outFilterMultimapNmax 1 \
	--outFilterMultimapScoreRange 2 \
	--outFilterMatchNminOverLread 0.95 \
	--outSAMmapqUnique 60 \

# outSAMmapqUnique is for being compatible with GATK

