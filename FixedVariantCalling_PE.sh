#!/bin/bash
#SBATCH --mem=144g
#SBATCH -n 24
#SBATCH -t 72:00:00
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=youremail@domain # send-to address
#SBATCH -o pipe.out
module load gcc/11.2.0
module load python
module load kraken
module load multiqc
module load fastp
module load fastqc
module load bowtie2
module load samtools
module load spades
module load singularity

# ---------- ADJUST THESE PARAMETERS ---------- #
lofreq=/path/to/lofreq
SRA_REPO=/path/to/srafiles
SNP_DB=/path/to/snpEffdb
WORK_DIR=/my/work/dir
THREADS=24
# ---------- ----------------------- ---------- #

TEMP_DIR=$WORK_DIR'/TEMP'
QC_DIR=$WORK_DIR'/QC'
FASTQ_DIR=$WORK_DIR/'FASTQs'
KRAKEN_DIR=$WORK_DIR/'KRAKEN'
VCF_DIR=$WORK_DIR'/VCFs'
CONS_DIR=$WORK_DIR'/Consensus_Sequences'
ASM_DIR=$WORK_DIR'/ASMs'
ASM_ALIGNMENTS=$ASM_DIR'/ALIGNMENTS'
ASM_INDICES=$ASM_DIR'/INDICES'
ASM_ANNOTATION=$ASM_DIR'/ANNOTATIONS'
ASM_BT2=$ASM_DIR'/BT2'

mkdir -p $TEMP_DIR
mkdir -p $QC_DIR
mkdir -p $FASTQ_DIR
mkdir -p $KRAKEN_DIR
mkdir -p $ASM_DIR
mkdir -p $ASM_ALIGNMENTS
mkdir -p $ASM_INDICES
mkdir -p $ASM_ANNOTATION
mkdir -p $ASM_BT2
mkdir -p $VCF_DIR
mkdir -p $CONS_DIR


#Download from SRA
#fasterq-dump SRR26891777 -e 12
### WORKFLOW ONE: GENERATE VD1 CONSENSUS GENOMES
while IFS=\" read -r line; do
if [[ $line == *"SRR"* ]]; then

    #STEP 1 - Build Sample Parameter Objects
    #"6","SRR","TX","DOSE","SN","GROUP"
    #"8","SRR","TX","DOSE","SN","GROUP"
    #"9","SRR","TX","DOSE","SN","GROUP"
    newl1=$(echo $line | awk -F\" '{print $4}') #SRR
    newl2=$(echo $line | awk -F\" '{print $6}') #Tx
    newl3=$(echo $line | awk -F\" '{print $8}') #Dose
    newl4=$(echo $line | awk -F\" '{print $10}') #PID
    newl5=$(echo $line | awk -F\" '{print $12}') #VisitDay
    out=$newl1"_"$newl2"_"$newl3"_"$newl4"_"$newl5
    ASM_NAME=$newl2"_"$newl3"_"$newl4

    #STEP 2 - QC reads, base quality sliding window 5' + base correction overlap + Overrepresented Seqs + Trim Poly X tail
    for r1 in $SRA_REPO'/'$newl1'_1.fastq'
    do
        r2=$(echo $r1 | sed -e 's/_1.fastq/_2.fastq/g')
        fastp -5 -c --overlap_diff_limit 0 -p --trim_poly_x -w $THREADS --length_required 120 -i $r1 -I $r2 \
            -o $TEMP_DIR'/'$newl1"_1.trim.fastq" -O $TEMP_DIR'/'$newl1"_2.trim.fastq" -j $QC_DIR'/'$newl1".fastp.json" -h $QC_DIR'/'$newl1".fastp.html"
        sed -E 's/([@+][A-Z]+[0-9]+\.[0-9]+) /\1\/1 /' $TEMP_DIR'/'$newl1"_1.trim.fastq" > $TEMP_DIR'/'$newl1"_1.trim.mod.fastq"
        sed -E 's/([@+][A-Z]+[0-9]+\.[0-9]+) /\1\/2 /' $TEMP_DIR'/'$newl1"_2.trim.fastq" > $TEMP_DIR'/'$newl1"_2.trim.mod.fastq"
    done
    fastqc -t $THREADS $TEMP_DIR'/'$newl1"_1.trim.mod.fastq" -o $QC_DIR
    fastqc -t $THREADS $TEMP_DIR'/'$newl1"_2.trim.mod.fastq" -o $QC_DIR
    multiqc -o $QC_DIR'/MQC' $QC_DIR

    #STEP 3 - Classify reads and filter metagenome to reads classified as Sarbecovirus -t=2509511
    var1=$TEMP_DIR'/'$newl1'_1.trim.mod.fastq'
    var2=$TEMP_DIR'/'$newl1'_2.trim.mod.fastq'
    kraken2 --use-names --paired --threads $THREADS --db /path/to/kraken2/2.1.2/KRAKEN2-HS-FUNGI-MM-DB \
        --report $KRAKEN_DIR'/'$newl1".report" $var1 $var2 > $KRAKEN_DIR'/'$newl1".kraken"
    python3 /path/to/krakenextract.py -k $KRAKEN_DIR'/'$newl1".kraken" \
        -s $var1 -s2 $var2 --report $KRAKEN_DIR'/'$newl1".report" \
        -t 2509511 --include-children --fastq-output \
        -o $FASTQ_DIR'/'$newl1'_1_filtered.fastq' -o2 $FASTQ_DIR'/'$newl1'_2_filtered.fastq'

    #STEP 4 - Align reads to reference | SARS2 for VD1 reference
    bowtie2 --phred33 -q -p $THREADS -x /sars2/bt2/reference/ \
        -1 $FASTQ_DIR'/'$newl1'_1_filtered.fastq' -2 $FASTQ_DIR'/'$newl1'_2_filtered.fastq' -S $ASM_ALIGNMENTS'/'$out".output.sam"
    samtools view -S -b $ASM_ALIGNMENTS'/'$out".output.sam" > $ASM_ALIGNMENTS'/'$out".output.bam"
    samtools sort $ASM_ALIGNMENTS'/'$out".output.bam" > $ASM_ALIGNMENTS'/'$out".sorted.bam"

    #STEP 5 - Assemble viral reads using SPADES with option --corona
    spades.py --rnaviral --corona -1 $FASTQ_DIR'/'$newl1'_1_filtered.fastq' -2 $FASTQ_DIR'/'$newl1'_2_filtered.fastq' \
        -o $ASM_DIR'/'$ASM_NAME"_spades_assembly" -t $THREADS -m 144
    cp ~/python_lib/cut_scaffolds.py $ASM_DIR
    python3 $ASM_DIR'/cut_scaffolds.py'

    #STEP 6 - Prokka Annotation + Bowtie2 Index
    singularity exec ~/prokka.sif prokka --kingdom Viruses \
        --outdir $ASM_ANNOTATION'/'$ASM_NAME --prefix $ASM_NAME --force --cpus $THREADS \
        --proteins ~/tryseqs.fasta $ASM_DIR'/'$ASM_NAME"_spades_assembly/scaffolds.mod.fasta"
    bowtie2 build -p $THREADS $ASM_ANNOTATION'/'$ASM_NAME'/'$ASM_NAME".fna" $ASM_BT2'/'$ASM_NAME

    #STEP 7 - Build snpEff database from prokka annotated GFFs
    toappend1='#'$ASM_NAME
    toappend2='\n'$ASM_NAME'.genome : '$ASM_NAME
    echo $toappend1 >> $SNP_DB'/snpEff.config'
    echo $toappend2 >> $SNP_DB'/snpEff.config'
    mkdir $SNP_DB'/GFFs/'$ASM_NAME
    cp $ASM_ANNOTATION'/'$ASM_NAME'/'$ASM_NAME".gff" $SNP_DB'/GFFs/'$ASM_NAME/'genes.gff'
    snpEff build -c $SNP_DB'/snpEff.config' -gff3 -v $ASM_NAME
fi
done <<< $(cat redux.csv)

### WORKFLOW 2: VARIANT CALLING OF POST VD1 SAMPLES

while IFS=\" read -r line; do
if [[ $line == *"SRR"* ]]; then

    #STEP 1 - Build Sample Parameter Objects
    newl1=$(echo $line | awk -F\" '{print $4}') #SRR
    newl2=$(echo $line | awk -F\" '{print $6}') #Tx
    newl3=$(echo $line | awk -F\" '{print $8}') #Dose
    newl4=$(echo $line | awk -F\" '{print $10}') #PID
    newl5=$(echo $line | awk -F\" '{print $12}') #VisitDay
    out=$newl1"_"$newl2"_"$newl3"_"$newl4"_"$newl5
    ASM_NAME=$newl2"_"$newl3"_"$newl4

    #STEP 2 - QC reads, base quality sliding window 5' + base correction overlap + Overrepresented Seqs + Trim Poly X tail
    for r1 in $SRA_REPO'/'$newl1'_1.fastq'
    do
        r2=$(echo $r1 | sed -e 's/_1.fastq/_2.fastq/g')
        fastp -5 -c --overlap_diff_limit 0 -p --trim_poly_x -w $THREADS --length_required 120 -i $r1 -I $r2 \
            -o $TEMP_DIR'/'$newl1"_1.trim.fastq" -O $TEMP_DIR'/'$newl1"_2.trim.fastq" -j $QC_DIR'/'$newl1".fastp.json" -h $QC_DIR'/'$newl1".fastp.html"
        sed -E 's/([@+][A-Z]+[0-9]+\.[0-9]+) /\1\/1 /' $TEMP_DIR'/'$newl1"_1.trim.fastq" > $TEMP_DIR'/'$newl1"_1.trim.mod.fastq"
        sed -E 's/([@+][A-Z]+[0-9]+\.[0-9]+) /\1\/2 /' $TEMP_DIR'/'$newl1"_2.trim.fastq" > $TEMP_DIR'/'$newl1"_2.trim.mod.fastq"
    done
    fastqc -t $THREADS $TEMP_DIR'/'$newl1"_1.trim.mod.fastq" -o $QC_DIR
    fastqc -t $THREADS $TEMP_DIR'/'$newl1"_2.trim.mod.fastq" -o $QC_DIR
    multiqc -o $QC_DIR'/MQC' $QC_DIR

    #STEP 3 - Classify reads and filter metagenome to reads classified as Sarbecovirus -t=2509511
    var1=$TEMP_DIR'/'$newl1'_1.trim.mod.fastq'
    var2=$TEMP_DIR'/'$newl1'_2.trim.mod.fastq'
    kraken2 --use-names --paired --threads $THREADS --db /path/to/kraken2/2.1.2/KRAKEN2-HS-FUNGI-MM-DB \
        --report $KRAKEN_DIR'/'$newl1".report" $var1 $var2 > $KRAKEN_DIR'/'$newl1".kraken"
    python3 /path/to/krakenextract.py -k $KRAKEN_DIR'/'$newl1".kraken" \
        -s $var1 -s2 $var2 --report $KRAKEN_DIR'/'$newl1".report" \
        -t 2509511 --include-children --fastq-output \
        -o $FASTQ_DIR'/'$newl1'_1_filtered.fastq' -o2 $FASTQ_DIR'/'$newl1'_2_filtered.fastq'

    #STEP 4 - Align reads to reference | SARS2 for VD1 reference
    bowtie2 --phred33 -q -p $THREADS -x $ASM_BT2'/'$ASM_NAME \
        -1 $FASTQ_DIR'/'$newl1'_1_filtered.fastq' -2 $FASTQ_DIR'/'$newl1'_2_filtered.fastq' -S $ASM_ALIGNMENTS'/'$out".output.sam"
    samtools view -S -b $ASM_ALIGNMENTS'/'$out".output.sam" > $ASM_ALIGNMENTS'/'$out".output.bam"
    samtools sort $ASM_ALIGNMENTS'/'$out".output.bam" > $ASM_ALIGNMENTS'/'$out".sorted.bam"

    #STEP 5 - Variant Calling Preprocessing
    $lofreq indelqual --dindel --ref $ASM_ANNOTATION'/'$ASM_NAME'/'$ASM_NAME".fna" -o $ASM_ALIGNMENTS'/'$out'.lofreq.bam' $ASM_ALIGNMENTS'/'$out'.sorted.bam'
    samtools index -o $ASM_ALIGNMENTS'/'$out'.lofreq.bam.bai' $ASM_ALIGNMENTS'/'$out'.lofreq.bam'

    #STEP 5.5 Variant Calling Process + Filter
    $lofreq call-parallel --call-indels --no-default-filter --pp-threads $THREADS -f $ASM_ANNOTATION'/'$ASM_NAME'/'$ASM_NAME".fna" --min-cov 10 \
        -o $ASM_ALIGNMENTS'/'$out'.lofreq.raw.vcf' $ASM_ALIGNMENTS'/'$out'.lofreq.bam'
    $lofreq filter --no-defaults --af-min 0.05 -i $ASM_ALIGNMENTS'/'$out'.lofreq.raw.vcf' -o $VCF_DIR'/'$out'.lofreq.5af.filtered.vcf'
    bgzip $VCF_DIR'/'$out'.lofreq.5af.filtered.vcf'
    bcftools index --threads $THREADS $VCF_DIR'/'$out'.lofreq.5af.filtered.vcf.gz'

    #PROBABLY NOT NEEDED, SAFE TO IGNORE /// STEP 6 Mask Low Coverage
    #bedtools genomecov -ibam $ASM_ALIGNMENTS'/'$out'.lofreq.bam' -bga | awk '$4 < 10 {print $0}' | \
    #    awk '{{print($1 "\t" $2 + 1 "\t" $3 "\tlow_coverage")}}' |\
    #    bedtools subtract -a - -b $VCF_DIR'/'$out'.lofreq.5af.filtered.vcf.gz' > $ASM_ALIGNMENTS'/'$out'_noncov.bed'

    #PROBABLY NOT NEEDED, SAFE TO IGNORE /// STEP 7 - Consensus Sequence
    #bcftools consensus -f $ASM_ANNOTATION'/'$ASM_NAME'/'$ASM_NAME".fna" \
    #-m $ASM_ALIGNMENTS'/'$out'_noncov.bed' -s - $VCF_DIR'/'$out'.lofreq.5af.filtered.vcf.gz' > $CONS_DIR'/'$out"_consensus.fasta"

    #STEP 6 Annotate VCFs with snpEff
    snpEff $ASM_NAME $VCF_DIR'/'$out'.lofreq.5af.filtered.vcf' > $VCF_DIR'/'$out'.FINAL_ANNOTATED.vcf'
fi
done <<< $(cat VD[X].csv)

