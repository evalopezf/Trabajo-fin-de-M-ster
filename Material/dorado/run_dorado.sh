outdir="/pool_sas/crh/analysis/twins/basecall_SUP"

CUDA_VISIBLE_DEVICES=0 /home/evalo/dorado-0.7.0-linux-x64/bin/dorado basecaller \

/home/evalo/models_dorado/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \

/pool_sas/crh/raw_data/twins/21-061_pod5/2023-09_GEMELOS_21-061.pod5 \

--reference /home/evalo/References/hg38.fa \

--modified-bases-models /home/evalo/models_dorado/dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mCG_5hmCG@v3.1 \

--secondary no \

--min-qscore 10 > $outdir/2023-09_GEMELOS_21-061.bam 2>$outdir/2023-09_GEMELOS_21-061_error.txt

if [ $? -eq 0 ]; then

samtools sort $outdir/2023-09_GEMELOS_21-061.bam -o $outdir/2023-09_GEMELOS_21-061_sorted.bam

samtools index $outdir/2023-09_GEMELOS_21-061_sorted.bam

else

echo "dorado basecaller failed for sample 2023-09_GEMELOS_21-061, check $outdir/2023-09_GEMELOS_21-061_error.txt for details"

fi
