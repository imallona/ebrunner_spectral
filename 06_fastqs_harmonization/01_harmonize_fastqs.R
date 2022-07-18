#!/usr/bin/env/R
##
## Preprocessing/harmonizing of bdrhapsody beads prepended seqs
##   Requires a R1 fastq (gz compressed) input (see help)
## Outputs a folder with three files
##   - non_matching.fastq.gz, reads that don't have the expected linkers/structure
##   - GTGA_GACA_linkers.fastq.gz, the oligodT probes
##   - AATG_CCAC_linkers.fastq.gz, the TSO probes
##
## Usage:
## $ Rscript 01_harmonize_fastqs.R --help
## usage: 01_harmonize_fastqs.R [-h] [-r1 r1] [-r2 r2] [-o output]

## optional arguments:
##   -h, --help            show this help message and exit
##   -r1 r1, --read1 r1    Path to the R1 fastq file (gz compressed); this file
##                         contains the barcodes
##   -r2 r2, --read2 r2    Path to the R2 fastq file (gz compressed); this file
##                         contains the cDNA
##   -o output, --output output
##                         Output folder name (creates the folder and overwrites
##                         content)

##  rationale: reads are a mixture of:
##  empty - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - T 25nt
##  empty - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - TATGCGTAGTAGGTATG
##  A     - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - T 25nt
##  A     - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - TATGCGTAGTAGGTATG
##  GT    - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - T 25nt
##  GT    - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - TATGCGTAGTAGGTATG
##  TCA   - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - T 25nt
##  TCGA  - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - TATGCGTAGTAGGTATG
##
## So we match a regex to capture linker1 and linker2, and remove the (potentially) prepended first 0-3 nts
##
## Izaskun Mallona izaskun.mallona@gmail.com
## GPLv3
##
## 04 July 2022

## .libPaths('/home/imallona/R/R4_bioc314')

suppressPackageStartupMessages({
    library("argparse")
    library('R.utils')
})

parser <- ArgumentParser()

parser$add_argument("-r1", "--read1", type="character", default='none',
                    help="Path to the R1 fastq file (gz compressed); this file contains the barcodes",
                    metavar="r1")

parser$add_argument("-r2", "--read2", type="character", default='none',
                    help="Path to the R2 fastq file (gz compressed); this file contains the cDNA",
                    metavar="r2")

parser$add_argument("-o", "--output", type="character", default='output',
                    help="Output folder name (creates the folder and overwrites content)",
                    metavar="output")

tryCatch({
    args <- parser$parse_args()
}, error = function(x) parser$print_help())


if (args$read1 != 'none' | args$read2 != 'none') {
    cat(sprintf('Processing started \t%s\n   R1: %s\n   R2: %s\n   Output folder: %s\n',
                Sys.time(),
                args$read1, args$read2, args$output))
} else {
    parser$print_help()
    stop('Malformed command, fastqs missing')
}


dir.create(args$output, showWarnings = FALSE, recursive = TRUE)

regex <- '([ACTGN]{0,3})([ACTGN]{9})(GTGA|AATG{1})([ACTGN]{9})(GACA|CCAC){1}([ACTGN]{9})([ACTGN]{8})(.*)'

## this is just for finetuning the regex; commented out
if (FALSE) {
    test1 <- c('TNAATGTAATGGGTGAACGACCACCGACAAAGGGAACTTCTTGATGTTTTTTTTTTTTTT',
               'CNCATTGCAAATGAATCCTGAACCACTTCAGCTCAGTTAAATATATGCGTAGGAGGTATG',
               'GNACACACAAAGTGAAGATAGTTCGACAACACCTTAGCTGCTTACTTTTTTTGTTTTTTT',
               'ANTTAGATGAATGCAGAAATCGCCACGATGGTCCACAGTTTCTTATGCCTACAACCGATG')

    test2 <- c('@A01251:431:H33CTDRX2:1:2101:30264:1016 1:N:0:GCTACGCT',
               'TNATAGCTTGTAGTGAAGGTTCGCTGACATGCATAGTAGACAGCGGTTTTTTTTTTTTTT',
               '+',
               'F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF:F,FFFFF')

    x <- strcapture(pattern = regex,
                    x = test1,
                    perl = TRUE,
                    proto = list(prepend = character(),
                                 seq1 = character(),
                                 link1 = character(),
                                 seq2 = character(),
                                 link2 = character(),
                                 seq3 = character(),
                                 umi = character(),
                                 tail = character()))
}


## extracts components of a fastq stanza and removes the prepended nts
## returns a tokenized fastq stanza, looking like this
##   prepend      seq1 link1      seq2 link2      seq3      umi              tail
## 1     TNA ATGTAATGG  GTGA ACGACCACC  GACA AAGGGAACT TCTTGATG    TTTTTTTTTTTTTT
tokenize_stanza <- function(original_stanza, regex) {
    x <- strcapture(pattern = regex,
            x = original_stanza[2],
            perl = TRUE,
            proto = list(prepend = character(),
                         seq1 = character(),
                         link1 = character(),
                         seq2 = character(),
                         link2 = character(),
                         seq3 = character(),
                         umi = character(),
                         tail = character()))

    x$seq_trimmed <- paste(c(x$seq1, x$link1, x$seq2, x$link2, x$seq3, x$umi, x$tail), collapse = '')
    x$qual_trimmed <- substr(original_stanza[4], nchar(x$prepend) + 1, nchar(original_stanza[4]))

    return(x)
}


## original stanza is the original fastq, a string `\n` separated
## x is the result of capturing the regex groups, as generated by tokenize_stanza
reconstruct_stanza <- function(original_stanza, x) {
    trimmed_stanza <- sprintf('%s\n%s\n%s\n%s',
                              original_stanza[1],
                              x$seq_trimmed,
                              original_stanza[3],
                              x$qual_trimmed)
    return(trimmed_stanza)
}

## being x a tokenized result, as generated by tokenize_stanza
count_number_ts_in_tail <- function(x) {
    return(sum(unlist(strsplit(x$tail, split = '')) == 'T'))
}

## being x a tokenized result, as generated by tokenize_stanza
## the number of Ts here is arbitrary and harcoded (5 consecutive Ts)
has_TTTTT_in_tail <- function(x) {
      return(grepl(x = x$tail, pattern = 'TTTTT'))
}


## being x a tokenized result, as generated by tokenize_stanza
## returns a boolean: has a TSO or no
## it allows 1 subst, 1 insertion, and as many deletions as the length of the prepended token
#    because, often, TSO is not fully seq in its 3'
has_tso_in_tail <- function(x) {
    tso <- 'TATGCGTAGTAGGTATG'

    ## x <- list(tail = 'TATGCGTAGTAGGTGA',
    ##           prepend = 'AT')

    ## x <- list(tail = 'TATGCGTAGTATGGGTGA',
    ##           prepend = '')

    ## agrep(x = tso, pattern = x$tail, max.distance = list(substitutions = 1,
    ##                                                             insertions = 1,
    ##                                                      deletions = nchar(x$prepend)))

    ## @todo explore further what's x, and what's pattern; better the way around?
    return(agrepl(x = tso, pattern = x$tail, max.distance = list(substitutions = 1,
                                                                 insertions = 1,
                                                                 deletions = nchar(x$prepend))))
}


initialize_counts_object <- function() {
    ## a log/summary to do some descriptive stats
    counts <- list('non_matching' = 0,
                   'GTGA' = list('GACA' = list('TSO' = 0,
                                               'Tn' = as.list(setNames(rep(0, 25), sprintf('T%s', 1:25)))),
                                 'else' = 0),
                   'AATG' = list('CCAC' = list('TSO' = 0,
                                               'Tn' = as.list(setNames(rep(0, 25), sprintf('T%s', 1:25)))),
                                 'else' = 0))
    return(counts)
}

## original_stanza a barcode fastq stanza (4 lines, `\n` separeted)
## original-cdna, the matching cDNA fastq stanza
## x is the tokenized result after capturing the regex groups
## output_dir the output directory path
## counts, an object with some stats
sort_reads_by_regex_match <- function(original_stanza, original_cdna, x, output_dir, counts) {

    
    ## Split reads into two buckets: TSO and not TSO.
    ## We do that according to the link1 and link2 sequences, and not checking whether the tail looks good or not
    ##   (even though TSO should have TSO seqs, and non TSO should have a polymer of Ts)

    ## the read does not match the regex
    if (any(is.na(x))) {
        write(original_stanza, file = file.path(output_dir, 'non_matching_R1.fastq'), append = TRUE)
        write(original_cdna, file = file.path(output_dir, 'non_matching_R2.fastq'), append = TRUE)

        counts[['non_matching']] <- counts[['non_matching']] + 1

    }
    else {
        ## read is supposed to belong to bucket oligodT, but we check the tail
        if (x$link1 == 'GTGA' & x$link2 == 'GACA') {
            if (has_tso_in_tail(x = x)) {
               
                write(reconstruct_stanza(x = x, original_stanza = original_stanza),
                      file = file.path(output_dir, 'GTGA_GACA_TSO_R1.fastq'), append = TRUE)

                write(original_cdna,
                      file = file.path(output_dir, 'GTGA_GACA_TSO_R2.fastq'), append = TRUE)

                counts[['GTGA']][['GACA']][['TSO']] <- counts[['GTGA']][['GACA']][['TSO']] + 1
            }
            else if (has_TTTTT_in_tail(x)) {
                num_ts <- sprintf('T%s', count_number_ts_in_tail(x))

                write(reconstruct_stanza(x = x, original_stanza = original_stanza),
                      file = file.path(output_dir, 'GTGA_GACA_Tn_R1.fastq'), append = TRUE)

                write(original_cdna,
                      file = file.path(output_dir, 'GTGA_GACA_Tn_R2.fastq'), append = TRUE)

                counts[['GTGA']][['GACA']][['Tn']][[num_ts]] <-counts[['GTGA']][['GACA']][['Tn']][[num_ts]] + 1
                
            }
            else {
                num_ts <- sprintf('T%s', count_number_ts_in_tail(x))

                write(reconstruct_stanza(x = x, original_stanza = original_stanza),
                      file = file.path(output_dir, 'GTGA_GACA_else_R1.fastq'), append = TRUE)

                write(original_cdna,
                      file = file.path(output_dir, 'GTGA_GACA_else_R2.fastq'), append = TRUE)

                ## we count the Ts anyway
                counts[['GTGA']][['GACA']][['Tn']][[num_ts]] <-counts[['GTGA']][['GACA']][['Tn']][[num_ts]] + 1
            }
        }
        else if (x$link1 == 'AATG' & x$link2 == 'CCAC') {
            if (has_tso_in_tail(x = x)) {
                
                write(reconstruct_stanza(x = x, original_stanza = original_stanza),
                      file = file.path(output_dir, 'AATG_CCAC_TSO_R1.fastq'), append = TRUE)

                write(original_cdna,
                      file = file.path(output_dir, 'AATG_CCAC_TSO_R2.fastq'), append = TRUE)

                counts[['AATG']][['CCAC']][['TSO']] <- counts[['AATG']][['CCAC']][['TSO']] + 1
            }
            else if (has_TTTTT_in_tail(x)) {
                num_ts <- sprintf('T%s', count_number_ts_in_tail(x))

                write(reconstruct_stanza(x = x, original_stanza = original_stanza),
                      file = file.path(output_dir, 'AATG_CCAC_Tn_R1.fastq'), append = TRUE)

                write(original_cdna,
                      file = file.path(output_dir, 'AATG_CCAC_Tn_R2.fastq'), append = TRUE)

                counts[['AATG']][['CCAC']][['Tn']][[num_ts]] <-counts[['AATG']][['CCAC']][['Tn']][[num_ts]] + 1
                
            }
            else {
                num_ts <- sprintf('T%s', count_number_ts_in_tail(x))

                write(reconstruct_stanza(x = x, original_stanza = original_stanza),
                      file = file.path(output_dir, 'AATG_CCAC_else_R1.fastq'), append = TRUE)

                write(original_cdna,
                      file = file.path(output_dir, 'AATG_CCAC_else_R2.fastq'), append = TRUE)

                counts[['AATG']][['CCAC']][['Tn']][[num_ts]] <-counts[['AATG']][['CCAC']][['Tn']][[num_ts]] + 1
            }

        }
        else{
            ## read matches the regex but the two linkers are not as expected (?)
            ## this use case might not be needed
            ## we won't count the Ts, neither
            write(original_stanza, file = file.path(output_dir, 'other_linkers_R1.fastq'), append = TRUE)
            write(original_cdna, file = file.path(output_dir, 'other_linkers_R2.fastq'), append = TRUE)
        }
    }
    return(counts)
}


## iterates over each stanza from a fastq gz connection
process_stanzas <- function(barcode_fn, cdna_fn, output_dir, regex) {
    ## the input filehandle
    fh <- gzfile(barcode_fn, 'r')
    cdna_fh <- gzfile(cdna_fn, 'r')

    ## the output files (barcodes, R1-like)
    nonmatching <- file(file.path(output_dir, 'non_matching_R1.fastq'), 'w')
    dt_proper <- file(file.path(output_dir, 'GTGA_GACA_Tn_R1.fastq'), 'w')
    tso_proper <- file(file.path(output_dir, 'AATG_CCAC_TSO_R1.fastq'), 'w')
    dt_not <- file(file.path(output_dir, 'GTGA_GACA_TSO_R1.fastq'), 'w')
    tso_not <- file(file.path(output_dir, 'AATG_CCAC_Tn_R1.fastq'), 'w')
    dt_else <- file(file.path(output_dir, 'GTGA_GACA_else_R1.fastq'), 'w')
    tso_else <- file(file.path(output_dir, 'AATG_CCAC_else_R1.fastq'), 'w')
    unusual <- file(file.path(output_dir, 'other_linkers_R1.fastq'), 'w')
        
    ## the output cDNA files (R2)
    nonmatching_cdna <- file(file.path(output_dir, 'non_matching_R1.fastq'), 'w')
    dt_proper_cdna <- file(file.path(output_dir, 'GTGA_GACA_Tn_R1.fastq'), 'w')
    tso_proper_cdna <- file(file.path(output_dir, 'AATG_CCAC_TSO_R1.fastq'), 'w')
    dt_not_cdna <- file(file.path(output_dir, 'GTGA_GACA_TSO_R1.fastq'), 'w')
    tso_not_cdna <- file(file.path(output_dir, 'AATG_CCAC_Tn_R1.fastq'), 'w')
    dt_else_cdna <- file(file.path(output_dir, 'GTGA_GACA_else_R1.fastq'), 'w')
    tso_else_cdna <- file(file.path(output_dir, 'AATG_CCAC_else_R1.fastq'), 'w')
    unusual_cdna <- file(file.path(output_dir, 'other_linkers_R1.fastq'), 'w')

    counts <- initialize_counts_object()
    
    while (TRUE) {
        original_barcode <- readLines(fh, n = 4)
        original_cdna <- readLines(cdna_fh, n = 4)
        
        if (length(original_barcode) == 0 | length(original_cdna) == 0) {
            break
        }
        
        counts <- sort_reads_by_regex_match(original_stanza = original_barcode,
                                  original_cdna = original_cdna,
                                  x = tokenize_stanza(original_stanza = original_barcode, regex = regex),
                                  output_dir = output_dir,
                                  counts = counts)
    }
    
    close(fh)
    close(cdna_fh)

    ## compress the fastq outputs
    for (item in c('nonmatching', 'dt_proper', 'tso_proper',
                   'dt_not', 'tso_not', 'dt_else', 'tso_else', 'unusual',
                   'nonmatching_cdna', 'dt_proper_cdna', 'tso_proper_cdna',
                   'dt_not_cdna', 'tso_not_cdna', 'dt_else_cdna', 'tso_else_cdna', 'unusual_cdna')) {
        close(get(item))
    }
    
    ## compress the fastq outputs
    for (item in list.files(output_dir, pattern = "*fastq")){
        ## close(file.path(output_dir, item))
        gzip(file.path(output_dir, item), ext = 'gz', overwrite = TRUE)
    }

    ## write the counts file
    saveRDS(object = counts, file = file.path(output_dir, 'counts.rds'))
}

process_stanzas(barcode_fn = args$read1,
                cdna_fn = args$read2, 
                output_dir = args$output,
                regex = regex)

cat(sprintf('End %s\t%s\n', args$output, Sys.time()))
