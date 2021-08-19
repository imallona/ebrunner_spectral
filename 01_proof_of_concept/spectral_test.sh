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

WD=/home/imallona/giulia                   # working dir, in portmac
GTF=gencode.vM27.basic.annotation.gff3.gz  # genes gtf
FEATURES=features.txt                      # genesymbols/ensg identifers of targets 
GENOME=GRCm39.primary_assembly.genome.fa   # genome gtf (not transcriptome)

# in case we don't have bedtools installed
export PATH=~/soft/bedtools/bedtools-2.29.2/bin/:"$PATH"

mkdir -p "$WD"; cd "$WD"

## get GTF
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/"$GTF"

## download genome
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/"$GENOME".gz

## get genes of interest, these are absolutely random
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

# get the 5' 300 bp ranges per gene (@todo check, is it per gene or per transcript?)
# for that, first get the chromosome sizes
gunzip "$GENOME".gz
samtools faidx "$GENOME"
cut -f1,2 "$GENOME".fai > mm39.chromsizes

## then, get these 300 bp after each gene start (strand aware of course)
grep -w gene selected.bed | \
    bedtools flank -i - -g mm39.chromsizes \
             -l 0 -r 300 -s > genes_300bp_up.bed

wc -l genes_300bp_up.bed

# remove these 5' regions from the regions to design the probes from
bedtools intersect -v -a  selected.bed \
         -b genes_300bp_up.bed | grep -w exon > safe_exons.bed

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

    bedtools getfasta -name -fi "$GENOME" -bed "$i".bed > "$i".fa
    wordcount --sequence "$i".fa -wordsize=25 -outfile "$i".wordcount

    # get only kmers appearing once within that exon
    grep -P "\t1$" "$i".wordcount | \
        awk -v id="$identifier" '{OFS=FS="\t"; print $0,id}' > \
            gene_"$gene"_trans_"$transcript"_exon_"$exon".out

    mkdir -p "$gene"
    mv gene_"$gene"_trans_"$transcript"_exon_"$exon".out "$gene"
    rm "$i".bed "$i".fa "$i".wordcount

    i=$(($i + 1))
done < safe_exons.bed

gzip "$GENOME"

# map (salmon/star/whatever)

# postprocess in R

## but with this strategy I don't know where these kmers were extracted from, just
##  from which exon did they come from
## anyway I could get this data while mapping, during the 'uniqueness' check step
## not sure this would be enough to then 'evenly space' them
