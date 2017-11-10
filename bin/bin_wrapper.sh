bin_size=$1
out_dir=$2
script_full_path=$(dirname "$0")

ref_fasta=$3
if [[ -z $ref_fasta ]]
then 
	ref_fasta=~/bioinfoSoftware/hg19/ucsc.hg19.fasta
fi

mkdir -p $out_dir

echo "running binner..."
python $script_full_path/binner.py $bin_size $out_dir/bins_1indx.tsv #!!! getting some duplicated bins

echo "getting 0-indx bins..."
awk -F "\t" 'BEGIN{OFS = "\t"}{print $1,$2-1,$3,$4}' $out_dir/bins_1indx.tsv > $out_dir/bins_0indx.tsv

echo "computing gc..."
echo $ref_fasta
bedtools nuc -fi $ref_fasta -bed $out_dir/bins_0indx.tsv | awk -F "\t" '{print $1,$2+1,$3,$4,$8+$9}' | grep -v "^#" > $out_dir/GC_bins.tsv

echo "FINISHED! results written to: $out_dir/GC_bins.tsv"


