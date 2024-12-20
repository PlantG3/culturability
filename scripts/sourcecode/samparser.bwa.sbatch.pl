#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my ($samfilecol, $mem, $time, $parserScript, $indir, $outdir, $help);
my ($filefeature, $threads, $outscriptOnly);
my $partition = "no";
my ($insert_min, $insert_max, $min_iden, $max_mismatch_perc,
    $keep_sam, $max_tail_perc, $gap, $min_score, $samtoolsModule);
my $result = &GetOptions("mem=s" => \$mem,
                         "time=s" => \$time,
			 "parserScript=s" => \$parserScript,
			 "samtoolsModule=s" => \$samtoolsModule,
			 "threads=i" => \$threads,
			 "indir=s" => \$indir,
			 "insert_min=i" => \$insert_min,
			 "insert_max=i" => \$insert_max,
			 "min_iden=i" => \$min_iden,
			 "max_mismatch_perc=i" => \$max_mismatch_perc,
			 "max_tail_perc=i" => \$max_tail_perc,
			 "gap=i" => \$gap,
			 "min_score=i" => \$min_score,
			 "filefeature=s" => \$filefeature,
			 "outdir=s" => \$outdir,
			 "partition=s" => \$partition,
			 "keepPARSEsam" => \$keep_sam,
			 "outscriptOnly" => \$outscriptOnly,
			 "help|h" => \$help
);

#$min_iden --mm $max_mismatch_perc 100 --tail $max_tail_perc 100 --gap $gap --mappingscore $min_score


#print help information if errors occur
if ($help){
	&errINF;
	exit;
}

$mem = "36G" if (!defined $mem);
$threads = 1 if (!defined $threads);
$time = "0-23:00:00" if (!defined $time);
$insert_min = 100 if (!defined $insert_min);
$insert_max = 1000 if (!defined $insert_max);
$min_iden = 50 if (!defined $min_iden);
$max_mismatch_perc = 6 if (!defined $max_mismatch_perc);
$max_tail_perc = 4 if (!defined $max_tail_perc);
$gap = 0 if (!defined $gap);
$min_score = 40 if (!defined $min_score);
$outdir = "." if (!defined $outdir);
$filefeature = ".sam" if (!defined $filefeature);
#$partition = "ksu-plantpath-liu3zhen.q,batch.q,killable.q" if (!defined $partition);


my $module_load;
if (defined $samtoolsModule) {
	$module_load = sprintf("module load %s", $samtoolsModule);
} else {
	print STDERR "no samtools module. --samtoolsModule is required.\n";
	exit 1;
}

open (IN, "ls \"$indir\" -1 |");
while (<IN>){
	chomp;
	my $samfile = $_;
	if ($samfile =~ $filefeature){
		my $sample = $samfile;
		$sample =~ s/$filefeature//g;
		print "$samfile\n";
		print "$sample\n";
		my $outfile = $sample.".sbatch";
		my $out = $sample.".parse.sam";
		my $out_log = $sample.".parse.log";
		open(OUT, ">", $outfile) || die;
		print OUT "#!/bin/bash -l\n";
		print OUT "#SBATCH --mem-per-cpu=$mem\n";
		print OUT "#SBATCH --time=$time\n";
		print OUT "#SBATCH --nodes=1\n";
		print OUT "#SBATCH --ntasks-per-node=$threads\n";
		if ($partition ne "no") {
			print OUT "#SBATCH --partition=$partition\n";
		}
		print OUT "perl $parserScript -i $indir\/$samfile --insert $insert_min $insert_max \\";
		print OUT "--identical $min_iden --mm $max_mismatch_perc 100 --tail $max_tail_perc 100 \\";
		print OUT "--gap $gap --mappingscore $min_score 1\>$outdir\/$out 2\>$outdir\/$out_log\n";
		if (defined $samtoolsModule) {
			print OUT "$module_load\n";
		}
		print OUT "samtools view -bS $out -@ $threads -o $sample.tmp.bam\n";
		print OUT "samtools sort -@ $threads -o $sample.bam $sample.tmp.bam\n";
		print OUT "samtools index -@ $threads $sample.bam\n";
		print OUT "rm $sample.tmp.bam\n";
		if (!$keep_sam) {
			print OUT "rm $outdir\/$out\n";
		}
		close OUT;
		my $sbatch_cmd = sprintf("sbatch %s", $outfile);
		print "$sbatch_cmd\n";
		if ($outscriptOnly) {
			print "$sbatch_cmd\n";
		} else {
			system($sbatch_cmd);
		}
		
	}
} 
close IN;

sub errINF {
	# to be added 
}
