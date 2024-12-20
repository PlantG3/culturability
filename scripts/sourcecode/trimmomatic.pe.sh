#!/bin/bash
#SBATCH --jobname=$jobname
#SBATCH --output=$trim_out_folder/1o-%j_trim.out
#SBATCH --error=$trim_out_folder/1o-%j_trim.err
#SBATCH --cpus-per-task=$trim_thread_num
#SBATCH --mems-per-cpu=$trim_mem_per_cpu
#SBATCH --time=$trim_time

# trimming adaaptor and no quality trimming
# designed for GBS 1-step trimming

# pre-trimming
module load Java/1.8.0_192
current_version=$(echo $trimmomatic | sed 's/.*\///g' | sed 's/.jar//g')

# trimming
#each_sample=$(echo $pair1_file | sed "s/$pair1_feature//g")
#pair2_file=$(echo $pair1_file | sed "s/$pair1_feature/$pair2_feature/g")
	
out_base=$trim_out_folder"/"$prefix
out_pair1_file=$out_base."R1.pair.fq"
out_single1_file=$out_base."R1.single.fq"
out_pair2_file=$out_base."R2.pair.fq"
out_single2_file=$out_base."R2.single.fq"

	
echo "Trimmomatic version: "$current_version
daytime=`date +%m/%d/%y-%H:%M:%S`
echo "Start trimming @"$daytime

### trimmomatic PE, refer the trimmomatic manual for the detail
java -jar $trimmomatic PE \
	-threads $trim_thread_num \
	$in_folder"/"$pair1_file $in_folder"/"$pair2_file \
	$out_pair1_file $out_single1_file \
	$out_pair2_file $out_single2_file \
	ILLUMINACLIP:$adaptor_file:3:20:10:1:true \
	LEADING:3 TRAILING:3 \
	SLIDINGWINDOW:4:0 \
	MINLEN:$trim_min_read_len \

### finished
daytime=`date +%m/%d/%y-%H:%M:%S`
echo "finish trimming @"$daytime

rm $out_single1_file
rm $out_single2_file

