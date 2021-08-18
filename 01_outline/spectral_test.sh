#!/bin/bash

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

WD=~/giulia
GTF=gencode.vM27.basic.annotation.gff3.gz
FEATURES=features.txt
# TARGETS="exon\|gene\|transcript"
GENOME=GRCm39.primary_assembly.genome.fa

export PATH=~/soft/bedtools/bedtools-2.29.2/bin/:"$PATH"

mkdir -p "$WD"; cd "$WD"

## get GTF
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/"$GTF"
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/"$GENOME".gz

## get genes of interest
cat << EOF > "$FEATURES"
Actb
Gm7341
Lactb2
Cyp20a1
Abca12
EOF

# grep genes of interest in gtf and transform to bed6

zcat "$GTF" | grep -v "^#" |\
    grep -f "$FEATURES" | grep -e "exon\|gene\|transcript" | \
    awk '{OFS=FS="\t"; print $1,$4,$5,$9,$3,$7}'  > selected.bed

wc -l selected.bed

# get the 5' 300 bp ranges per gene (@todo check) and not per gene
# first, get chromsizes
## @todo better use faidx from the fasta
# mysql --user=genome --host=genome-mysql.soe.ucsc.edu \
#       -B \
#       -A -D mm39 -N -e 'select chrom,size from chromInfo' | grep -v "_" > mm39.chromsizes

gunzip "$GENOME".gz
samtools faidx "$GENOME"
cut -f1,2 "$GENOME".fai > mm39.chromsizes


grep -w gene selected.bed | \
    bedtools flank -i - -g mm39.chromsizes \
             -l 0 -r 300 -s > genes_300bp_up.bed

wc -l genes_300bp_up.bed

# remove these 5' regions from the regions to design the probes from

bedtools intersect -v -a  selected.bed \
         -b genes_300bp_up.bed | grep -w exon > safe_exons.bed


# extract the fasta seq for each gene

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
    # mkdir -p "$gene"/"$transcript" ; cd "$_"
    
    bedtools getfasta -name -fi "$GENOME" -bed "$i".bed > "$i".fa
    wordcount --sequence "$i".fa -wordsize=25 -outfile "$i".wordcount

    # get only kmers appearing once within that exon
    grep -P "\t1$" "$i".wordcount | \
        awk -v id="$identifier" '{OFS=FS="\t"; print $0,id}' > \
            gene_"$gene"_trans_"$transcript"_exon_"$exon".out

    mkdir -p "$gene"
    mv gene_"$gene"_trans_"$transcript"_exon_"$exon".out "$gene"
    rm "$i".bed "$i".fa "$i".wordcount
    # cd ../..
    
    i=$(($i + 1))
done < safe_exons.bed

gzip "$GENOME"

# postprocess in R

## but with this strategy I don't know where these kmers were extracted from


# now,  word count all possible 25nt slices for each one of the 'exon' records,
#   keeping track of the gene name
#   so then we can aggregate them to maximize the number of exons, and
#    minimize the number of genes (in R)


# loop per fasta record and generate a tmp record, from which the kmers are generated

wordcount tembl:u68037 -wordsize=25, maybe use jellyfish


# extract the fasta sequence

# and group by which genes and exons belong to
