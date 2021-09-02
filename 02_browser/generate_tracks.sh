#!/bin/bash

# UCSC browseable tracks

WD=/home/imallona/giulia

# filtered exons
EXONS="$WD"/safe_exons.bed

MASK="$WD"/transcripts_300bp_up.bed

mkdir -p "$WD"/tracks

## track with probes

zcat probes/*tsv.gz | grep -v "^gene_id" | \
    awk '{ 
          OFS=FS="\t"; 
          print $7,$8,$9,$1,"0",$10
         }' | gzip -c > "$WD"/tracks/probes_track.bed.gz

echo 'browser position chr6:52142645-52304303
track name=spectral_probes description="Spectral PoC: probes (25nt, unique)" visibility=1' | gzip -c > header.gz

zcat header.gz "$WD"/tracks/probes_track.bed.gz | gzip -c > "$WD"/tracks/probes_track.bed.gz.tmp
mv -f "$WD"/tracks/probes_track.bed.gz.tmp "$WD"/tracks/probes_track.bed.gz

## tracks with 'target' exons

awk '{ 
       OFS=FS="\t"; 
       print $1,$2,$3,$4,"0",$6
      }' "$EXONS" | gzip -c > "$WD"/tracks/exons_track.bed.gz

echo 'browser position chr6:52142645-52304303
track name=spectral_exons description="Spectral PoC: target exons" color=0,0,255 visibility=1' | gzip -c > header.gz

zcat header.gz "$WD"/tracks/exons_track.bed.gz | gzip -c > "$WD"/tracks/exons_track.bed.gz.tmp
mv -f "$WD"/tracks/exons_track.bed.gz.tmp "$WD"/tracks/exons_track.bed.gz

## track with the mask (areas within the 300 bp of each tx, probes were not designed)

awk '{ 
       OFS=FS="\t"; 
       print $1,$2,$3,".","0",$6
      }' "$MASK" | gzip -c > "$WD"/tracks/mask_track.bed.gz

echo 'browser position chr6:52142645-52304303
track name=spectral_mask description="Spectral PoC: mask" color=255,0,0 visibility=1' | gzip -c > header.gz

zcat header.gz "$WD"/tracks/mask_track.bed.gz | gzip -c > "$WD"/tracks/mask_track.bed.gz.tmp
mv -f "$WD"/tracks/mask_track.bed.gz.tmp "$WD"/tracks/mask_track.bed.gz

rm header.gz

rsync -avt "$WD"/tracks/*gz imlspenticton:/var/www/imallona

echo 'http://imlspenticton.uzh.ch/imallona/mask_track.bed.gz
http://imlspenticton.uzh.ch/imallona/exons_track.bed.gz
http://imlspenticton.uzh.ch/imallona/probes_track.bed.gz' > "$WD"/tracks/spectral_poc.txt

rsync -avt "$WD"/tracks/spectral_poc.txt imlspenticton:/var/www/imallona

# http://genome.ucsc.edu/cgi-bin/hgTracks?genome=mm39&hgt.customText=http://imlspenticton.uzh.ch/imallona/spectral_poc.txt
