#!/bin/bash
#SBATCH --account=rrg-rbeiko
#SBATCH --job-name=Musquash_anlaysis  # Job name
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --cpus-per-task=32          # Number of CPU cores to use
#SBATCH --time=05:00:00             # Time limit (hh:mm:ss)
#SBATCH --mem=200G                   # Memory limit
#SBATCH --output=qiime_%j.log       # Standard output and error log
#SBATCH --error=qiime_%j.err        # Error log


export TMPDIR=/lustre09/project/6070434/mfares/PROJ-Musquash/tmp/

# Load Qiime2 module
module load qiime2/2024.5
module load blast+/2.14.1

target_dir="/lustre09/project/6070434/mfares/PROJ-Musquash/Analysis/data"
marker="COI"
year="25_B"
execute_block=true


if [ "$execute_block" = true ]; then  
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$target_dir/Musq-${year}-${marker}.manifest" \
  --output-path "$target_dir/Musq-${year}-${marker}-Reads.qza" \
  --input-format PairedEndFastqManifestPhred33



 
if [ "$marker" = "mifish" ]; then
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$target_dir/Musq-${year}-${marker}-Reads.qza" \
  --o-table "$target_dir/Musq-${year}-${marker}-feature-table.qza" \
  --o-representative-sequences "$target_dir/Musq-${year}-${marker}-rep-seq.qza" \
  --p-trunc-len-f 126 \
  --p-trunc-len-r 119 \
  --p-trim-left-f 5 \
  --p-trim-left-r 5 \
  --p-n-reads-learn 1000000 \
  --p-n-threads 32 \
  --o-denoising-stats "$target_dir/Musq-${year}-${marker}-denoise_stats.qza" \
  --verbose

elif [ "$marker" = "COI" ]; then
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$target_dir/Musq-${year}-${marker}-Reads.qza" \
  --o-table "$target_dir/Musq-${year}-${marker}-feature-table.qza" \
  --o-representative-sequences "$target_dir/Musq-${year}-${marker}-rep-seq.qza" \
  --p-trunc-len-f 200 \
  --p-trunc-len-r 180 \
  --p-trim-left-f 30 \
  --p-trim-left-r 30 \
  --p-n-reads-learn 1000000 \
  --p-n-threads 32 \
  --o-denoising-stats "$target_dir/Musq-${year}-${marker}-denoise_stats.qza" \
  --verbose
fi



qiime alignment mafft --i-sequences  "$target_dir/Musq-${year}-${marker}-rep-seq.qza" \
  --p-n-threads 32 \
  --o-alignment "$target_dir/Musq-${year}-${marker}-rep-seq-algn" \
  --verbose
  
qiime alignment mask --i-alignment "$target_dir/Musq-${year}-${marker}-rep-seq-algn.qza" \
  --o-masked-alignment "$target_dir/Musq-${year}-${marker}-rep-seq-algn-msk" \
  --verbose


qiime phylogeny fasttree --i-alignment "$target_dir/Musq-${year}-${marker}-rep-seq-algn-msk.qza" \
  --p-n-threads 32 \
  --o-tree "$target_dir/Musq-${year}-${marker}-unrooted-tree" \
  --verbose

qiime phylogeny midpoint-root --i-tree "$target_dir/Musq-${year}-${marker}-unrooted-tree.qza" \
  --o-rooted-tree  "$target_dir/Musq-${year}-${marker}-midrooted_tree" \
  --verbose
 



qiime tools export \
  --input-path "$target_dir/Musq-${year}-${marker}-rep-seq.qza" \
  --output-path "$target_dir/Musq-${year}-${marker}-rep-seq"

mv $target_dir/Musq-${year}-${marker}-rep-seq/dna-sequences.fasta $target_dir/Musq-${year}-${marker}-rep-seq.fasta
rm -r $target_dir/Musq-${year}-${marker}-rep-seq

#fi


export BLASTDB=/lustre09/project/6070434/shared/Databases/NCBI/nt_euk

 
## run with 10 top hits

blastn -db nt_euk \
-query $target_dir/Musq-${year}-${marker}-rep-seq.fasta \
-max_target_seqs 10 \
-out $target_dir/Musq-${year}-${marker}-Blast-nt_euk-90.tsv \
-evalue 1e-20 \
-perc_identity 90 \
-qcov_hsp_perc 90 \
-outfmt "6 qseqid sseqid pident qcovs evalue length staxids ssciname sblastname scomname" \
-num_threads 8



if [ "$marker" = "mifish" ]; then
blastdb_path="/lustre09/project/6070434/mfares/databases/12S/MIDORI/MIDORI2_UNIQ_NUC_GB266_srRNA_BLAST"
elif [ "$marker" = "COI" ]; then
blastdb_path="/lustre09/project/6070434/mfares/databases/COI/MIDORI/MIDORI2_UNIQ_NUC_GB265_CO1_BLAST"
fi



blastn \
-db "$blastdb_path" \
-query $target_dir/Musq-${year}-${marker}-rep-seq.fasta \
-max_target_seqs 10 \
-out $target_dir/Musq-${year}-${marker}-Blast-MIDORI-90.tsv \
-evalue 1e-20 \
-perc_identity 90 \
-qcov_hsp_perc 90 \
-outfmt "6 qseqid sseqid pident qcovs evalue length staxids ssciname sblastname scomname" \
-num_threads 8


fi

if [ "$marker" = "COI" ]; then

# Load java module 
module load java/21.0.1
rdp_dir="/lustre09/project/6070434/mfares/PROJ-Musquash/rdp_classifier"
 
java -Xmx128g -Xms8g -jar $rdp_dir/rdp_classifier_2.13/dist/classifier.jar -t $rdp_dir/mydata/rRNAClassifier.properties -o $target_dir/Musq-${year}-${marker}-RDP-results $target_dir/Musq-${year}-${marker}-rep-seq.fasta

fi
