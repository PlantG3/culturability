perl ../sourcecode/bwa.sbatch.pl \
        --mem 4G --time 0-23:00:00 \
        --bwa_shell $path_to_bwa \
	--indir ../1_debarcoding/all \
	--outdir . \
	--db $ref \
        --fq1feature .R1.pair.fq --fq2feature .R2.pair.fq \
        --threads 8

