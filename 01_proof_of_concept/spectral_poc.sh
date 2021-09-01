#!/bin/bash
##
## Erich Brunner's +  Giulia Moro's Spectral data processor (PoC) 
## 
## Software reqs: bedtools, STAR, EMBOSS, pigz
##
## Izaskun Mallona
##
## License GPLv3
##
## Started 18 Aug 2021

## Criteria 
# 1.Specific for gene (in this case transcription factors)
# 2.Non – overlapping capture sequences 
# 3.Recognising as many transcript variants as possible
# 4.~ 300 bp from 5’ end 
# 5.In CDS
# 6.Evenly spaced in gene
# 7.GC content upstream of  probe binding
# 8.Thermodynamic properties probe binding 

# ## Importance criteria

# A.In CDS
# B.Specificity 
# C.Thermodynamic properties determined by the program itself
# D.GC content upstream of probe
# E.As many transcript variants as possible 
# F.Evenly spaced in gene (+ non-overlapping) in as many CDSs as possible connected to E 
# G.~ 300 bp from 5’ end CDS

WD=/home/imallona/giulia                   # working dir, in portmac
GTF=gencode.vM27.basic.annotation.gff3     # genes gtf
## genesymbols/ensg identifers of targets
FEATURES=/home/imallona/src/ebrunner_spectral/01_proof_of_concept/data/tf_mouse_gmoro.txt 
GENOME=GRCm39.primary_assembly.genome.fa   # genome gtf (not transcriptome)
NTHREADS=5                                 # number of cores
KMER_LENGTH=25                             # kmer length
POSTPROC_RSCRIPT=/home/imallona/src/ebrunner_spectral/01_proof_of_concept/spectral_poc_postprocess.R

# binaries
STAR=~/soft/star/STAR-2.7.3a/bin/Linux_x86_64/STAR   
# in case we don't have bedtools installed
export PATH=~/soft/bedtools/bedtools-2.29.2/bin/:"$PATH"

mkdir -p "$WD"; cd "$WD"

## get GTF
if [ ! -f "$GTF".gz ]; then
    echo "Download GTF"
    wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/"$GTF".gz
fi

## download genome
if [ ! -f "$GENOME".gz ]; then
    echo "Download genome"
    wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/"$GENOME".gz
fi

## get genes of interest, these are absolutely random
# cat << EOF > "$FEATURES"
# Actb
# Gm7341
# Lactb2
# Cyp20a1
# Abca12
# Adgrg1
# Ndrg2
# Dmkn
# EOF


# Check how many genes of interest are there
echo "Designing probes for: `wc -l $FEATURES` genes"


# grep genes of interest in gtf and transform to bed6
# TODO brainstorm if exons are indeed to be fetched, and not CDS
pigz --decompress -p "$NTHREADS" --stdout --keep "$GTF".gz | \
    LC_ALL=C fgrep -v "^#" | \
    LC_ALL=C grep -w -f "$FEATURES" | \
    LC_ALL=C grep -e "exon\|gene\|transcript" | \
    awk '{OFS=FS="\t"; print $1,$4,$5,$9,$3,$7}'  > selected.bed

echo "Number of features: `wc -l selected.bed`"

# get the 5' 300 bp ranges per transcript (these will be removed from the analysis)
# for that, first get the chromosome sizes
pigz --decompress -p "$NTHREADS" "$GENOME".gz

if [ ! -f mm39.chromsizes ]; then
    samtools faidx "$GENOME"
    cut -f1,2 "$GENOME".fai > mm39.chromsizes
fi

## then, get these 300 bp after each transcript start (strand aware of course)
# LC_ALL=C grep -w gene selected.bed | \
#     bedtools flank -i - -g mm39.chromsizes \
#              -l 0 -r 300 -s > genes_300bp_up.bed
awk '$5== "transcript" {print $0}' selected.bed | \
        bedtools flank -i - -g mm39.chromsizes \
                 -l 0 -r 300 -s > transcripts_300bp_up.bed

echo "Number of 300 bp regions to avoid: `wc -l transcripts_300bp_up.bed`"

# remove these 5' regions from the regions to design the probes from
bedtools intersect -v -a  selected.bed \
         -b transcripts_300bp_up.bed | grep -w exon > safe_exons.bed

echo "Number of features (5' 300 bp excluded): `wc -l safe_exons.bed`"
# loop to extract unique (per exon) 25-mers to be stored with metadata, including:
#  -gene they belong to
#  -transcript they belong to
#  -exon they belong to
#
# e.g. output looks like
# $ ls -1 ENSMUSG00000025937.7/*out
# ENSMUSG00000025937.7/gene_ENSMUSG00000025937.7_trans_ENSMUST00000027071.7_exon_ENSMUST00000027071.7:1.out
# ENSMUSG00000025937.7/gene_ENSMUSG00000025937.7_trans_ENSMUST00000027071.7_exon_ENSMUST00000027071.7:2.out
# ENSMUSG00000025937.7/gene_ENSMUSG00000025937.7_trans_ENSMUST00000027071.7_exon_ENSMUST00000027071.7:3.out
# ENSMUSG00000025937.7/gene_ENSMUSG00000025937.7_trans_ENSMUST00000027071.7_exon_ENSMUST00000027071.7:4.out
# ENSMUSG00000025937.7/gene_ENSMUSG00000025937.7_trans_ENSMUST00000027071.7_exon_ENSMUST00000027071.7:5.out
# ENSMUSG00000025937.7/gene_ENSMUSG00000025937.7_trans_ENSMUST00000027071.7_exon_ENSMUST00000027071.7:6.out
# ENSMUSG00000025937.7/gene_ENSMUSG00000025937.7_trans_ENSMUST00000027071.7_exon_ENSMUST00000027071.7:7.out

# whereas one of these files looks like (two first k-mers)
# $ head -2 ENSMUSG00000025937.7/gene_ENSMUSG00000025937.7_trans_ENSMUST00000027071.7_exon_ENSMUST00000027071.7:1.out
# ATCTACCTCGGTGACCGCCGGCACC	1	exon_coord=chr1:13730553-13730770-;ID=exon:ENSMUST00000027071.7:1;Parent=ENSMUST00000027071.7;gene_id=ENSMUSG00000025937.7;transcript_id=ENSMUST00000027071.7;gene_type=protein_coding;gene_name=Lactb2;transcript_type=protein_coding;transcript_name=Lactb2-201;exon_number=1;exon_id=ENSMUSE00000154319.7;level=2;protein_id=ENSMUSP00000027071.6;transcript_support_level=1;mgi_id=MGI:2442551;tag=basic,appris_principal_1,CCDS;ccdsid=CCDS14823.1;havana_gene=OTTMUSG00000021599.3;havana_transcript=OTTMUST00000051280.2
# TATCTACCTCGGTGACCGCCGGCAC	1	exon_coord=chr1:13730553-13730770-;ID=exon:ENSMUST00000027071.7:1;Parent=ENSMUST00000027071.7;gene_id=ENSMUSG00000025937.7;transcript_id=ENSMUST00000027071.7;gene_type=protein_coding;gene_name=Lactb2;transcript_type=protein_coding;transcript_name=Lactb2-201;exon_number=1;exon_id=ENSMUSE00000154319.7;level=2;protein_id=ENSMUSP00000027071.6;transcript_support_level=1;mgi_id=MGI:2442551;tag=basic,appris_principal_1,CCDS;ccdsid=CCDS14823.1;havana_gene=OTTMUSG00000021599.3;havana_transcript=OTTMUST00000051280.2


i=0
while IFS= read -r line
do
    echo "$line" > "$i".bed
    ## adding the chr coords of the exon is a bit unnecessary, because
    ##   we have them in the GTF; anyway, adding them
    identifier=$(awk '{print "exon_coord="$1":"$2"-"$3$6";"$4}' "$i".bed)
    exon=$([[ $identifier =~ .*exon:(.*)\;Parent.* ]] && echo "${BASH_REMATCH[1]}")
    gene=$([[ $identifier =~ .*gene_id=(.*)\;transcript_id.* ]] && echo "${BASH_REMATCH[1]}")
    transcript=$([[ $identifier =~ .*transcript_id=(.*)\;gene_type.* ]] && echo "${BASH_REMATCH[1]}")
    id_short=$(awk -v exon="$exon" '{print "exon_coord="$1":"$2"-"$3$6";"exon}' "$i".bed)
    
    bedtools getfasta -name -fi "$GENOME" -bed "$i".bed > "$i".fa
    wordcount --sequence "$i".fa -wordsize="$KMER_LENGTH" -outfile "$i".wordcount &> /dev/null

    # get only kmers appearing once within that exon
    grep -P "\t1$" "$i".wordcount | \
        awk -v id="$id_short" '{OFS=FS="\t"; print $0,id}' > \
            gene_"$gene"_trans_"$transcript"_exon_"$exon".out

    awk -v id="$id_short" '{OFS=FS="\t";print ">"id"_kmer_"NR"\n"$1}' \
        gene_"$gene"_trans_"$transcript"_exon_"$exon".out > \
        gene_"$gene"_trans_"$transcript"_exon_"$exon".fa

    mkdir -p tmp/"$gene"
    gzip gene_"$gene"_trans_"$transcript"_exon_"$exon".fa -c > tmp/"$gene"/gene_"$gene"_trans_"$transcript"_exon_"$exon".fa.gz
    
    rm "$i".bed "$i".fa "$i".wordcount gene_"$gene"_trans_"$transcript"_exon_"$exon".out \
       gene_"$gene"_trans_"$transcript"_exon_"$exon".fa

    i=$(($i + 1))
done < safe_exons.bed

# merge all genes' kmers, for mapping

find ./tmp -name "gene*fa.gz" | xargs zcat | gzip -c > kmers_"$KMER_LENGTH".fa.gz
rm -rf ./tmp

# map (salmon/star/whatever)

# first attempt with a gDNA + GTF sort of approach (STAR)
#  maybe to be replaced with a pure transcriptome sort of mapping, afterwards
# note the index kmer length fitting to the kmer length (sjdbOverhang kmer length -1)

if [ ! -f  indices/$(basename "$GTF" .gff3)_kmer_"$KMER_LENGTH"/SA ]; then
    echo 'indexing'

    mkdir -p indices/$(basename "$GTF" .gff3)_kmer_"$KMER_LENGTH"

    pigz --decompress -p "$NTHREADS" "$GTF".gz

    "$STAR" --runThreadN "$NTHREADS" \
            --runMode genomeGenerate \
            --genomeDir indices/$(basename "$GTF" .gff3)_kmer_"$KMER_LENGTH" \
            --genomeFastaFiles "$GENOME" \
            --sjdbGTFfile "$GTF" \
            --sjdbOverhang $(($KMER_LENGTH - 1))

    pigz -p "$NTHREADS" "$GENOME"
    pigz -p "$NTHREADS" "$GTF"
fi

# map

mkdir -p mapping

# with nsorted output: multiple alignments of a read are adjacent
# outSAMmultNmax 1: report only one alignment for multimappers
"$STAR" --genomeDir indices/$(basename "$GTF" .gff3)_kmer_"$KMER_LENGTH" \
        --runThreadN "$NTHREADS" \
        --readFilesIn kmers_"$KMER_LENGTH".fa.gz \
        --outFileNamePrefix mapping/kmers_"$KMER_LENGTH" \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outSAMmultNmax 1 \
        --outSAMattributes NH HI NM MD AS nM \
        --readFilesCommand zcat

# awk parsing the bam file: items to be reported:
#   read_name
#   map_chr
#   map_star
#   map_end
#   map_strand
#   NH Number of reported alignments that contain the query in the current record <---
#   HI Query hit index
#   NM edit distance to the reference
#   MD string encoding mismatched and deleted reference bases
#   AS id the local alignment score (paired for paired-end reads).
#   nM number of mismatches per (paired) alignment

# we grep NH=1 to get unique mappers (as defined by STAR)
# exon_coord=chr5:142889256-142889696-;ENSMUST00000167721.8:5_kmer_45

samtools view mapping/kmers_"$KMER_LENGTH"Aligned.out.bam  | \
    grep -w 'NH:i:1' | \
    awk '{OFS=FS="\t"; print $1,$3,$4,$6,$12,$13,$14,$15,$16,$17,$10}' | \
    pigz -p $NTHREADS --stdout > mapping/kmers_"$KMER_LENGTH".uniques.gz


# postprocess in R, by chromosome

mkdir -p probes

## sequential version
# for chrom in $(seq 1 19) X Y
# do
#     chrom="chr"$chrom
#     echo $chrom

#     # run=$(printf 'rmarkdown::render("%s")' "$result" "$size" "$name" "$visits" "$inbound" "$outbound")
#     # R -e "rmarkdown::render('"P
#     Rscript "$POSTPROC_RSCRIPT" "$chrom" probes/spectral_poc_"$chrom".tsv
#     gunzip probes/spectral_poc_"$chrom".tsv
# done

## in batches of $NTHREADS chromosomes
N=$NTHREADS
(
    ## chrY does not have much on it/skipping
    for chrom in $(seq 1 19) X
    do
        chrom="chr""$chrom"

        ((i=i%N)); ((i++==0)) && wait
        Rscript "$POSTPROC_RSCRIPT" "$chrom" probes/spectral_poc_"$chrom".tsv &
    done
)

gzip probes/*tsv

# knit the 'final' report
