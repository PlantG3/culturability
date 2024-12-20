#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my ($fqfilecol, $mem, $time, $bwa, $indir, $db, $outdir, $help);
my ($minMapQual, $fq1feature, $fq2feature, $threads);
my $result = &GetOptions("mem=s" => \$mem,
                         "time=s" => \$time,
                         "bwa_shell=s" => \$bwa,
			 "indir=s" => \$indir,
			 "db=s" => \$db,
			 "outdir=s" => \$outdir,
			 "fq1feature=s" => \$fq1feature,
			 "fq2feature=s" => \$fq2feature,
			 "minmapqual=i" => \$minMapQual,
			 "threads=i" => \$threads,
			 "help|h" => \$help
);

# print help information if errors occur:
if ($help) {
	&errINF;
	exit;
}

$mem = "36G" if (!defined $mem);
$time = "48:00:00" if (!defined $time);
$fq1feature = ".R1.pair.fq" if (!defined $fq1feature);
$fq2feature = ".R2.pair.fq" if (!defined $fq2feature);
$outdir = "." if (!defined $outdir);
$minMapQual = 40 if (!defined $minMapQual);
$threads = 8 if (!defined $threads);

open (IN,"ls \"$indir\" -1 |");
while (<IN>) {
	chomp;
	my $fqfile = $_;
	if ($fqfile =~ $fq1feature) {
		my $sample = $fqfile;
		#$sample =~ s/.*\///g;
		$sample =~ s/$fq1feature//g;
		my $fq1 = $fqfile;
		my $fq2 = $fq1;
		$fq2 =~ s/$fq1feature/$fq2feature/g;
		print "$fqfile\n";
		print "$sample\n";
		my $outfile = $sample.".sbatch";
		my $out = $sample.".sam";
		
		my $infq1 = $indir."/".$fq1;
		if ($fq1 =~ /gz$/) {
			$infq1 = "<(gunzip -c ".$indir."/".$fq1.")";
		}	
		
		my $infq2 = $indir."/".$fq2; 
		if ($fq2 =~ /gz$/) {
			$infq2 = "<(gunzip -c ".$indir."/".$fq2.")";
		}

		open(OUT, ">", $outfile) || die;
		print OUT "#!/bin/bash -l\n";
		print OUT "#SBATCH --mem-per-cpu=$mem\n";
		print OUT "#SBATCH --time=$time\n";
		print OUT "#SBATCH --nodes=1\n";
		print OUT "#SBATCH --ntasks-per-node=$threads\n";
		print OUT "$bwa mem -t $threads -T $minMapQual -S -R \'\@RG\\tID\:$sample\\tSM\:$sample\' $db $infq1 $infq2 \> $outdir\/$out\n";
		close OUT;
		my $sbatch_cmd = sprintf("sbatch %s", $outfile);
		print "$sbatch_cmd\n";
		system($sbatch_cmd);
	}
}

close IN;

sub errINF {
	# to be added (SL)
}
