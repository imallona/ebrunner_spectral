#!/usr/bin/env/R
##
## Preprocessing/harmonizing of bdrhapsody beads prepended seqs
##   Requires a R1 fastq (gz compressed) input (see help)
## Outputs the, uncompressed, labelled and potentially trimmed  R1 reads to stdout
##   and R2 (cDNA) reads to stderr (!)
## Read names are annotated as
##   - else, the noncanonical reads
##   - GTGA_GACA_WTA, the oligodT/WTA probes
##   - AATG_CCAC_TSO, the TSO probes
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
## 04 July 2022, updated 4th Jan 2023

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

tryCatch({
    args <- parser$parse_args()
}, error = function(x) parser$print_help())


if (args$read1 == 'none' | args$read2 == 'none') {
    parser$print_help()
    stop('Malformed command, fastqs missing')
}

if (FALSE) {
    args <- list(read1 = './test/test_R1.fastq.gz')
                 
}


regex <- '([ACTGN]{0,3})([ACTGN]{9})(GTGA|AATG)([ACTGN]{9})(GACA|CCAC)([ACTGN]{9})([ACTGN]{8})(.*)'

## this is just for finetuning the regex; commented out
if (FALSE) {
    test1 <- c('TNAATGTAATGGGTGAACGACCACCGACAAAGGGAACTTCTTGATGTTTTTTTTTTTTTT',
               'CNCATTGCAAATGAATCCTGAACCACTTCAGCTCAGTTAAATATATGCGTAGGAGGTATG',
               'GNACACACAAAGTGTAGATAGTTCGACAACACCTTAGCTGCTTACGTGTGTTGTTTTTTT',
               'ANTTAGATGAATGCAGAAATCGCCACGATGGTCCACAGTTTCTTATGCCTACAACCGATG',
               'TGTAATCACAATGTGCGTATCACAACTGCATTGAAGACACTGCGGTGTGCTCCCAAAAAA')

    ## so the regex might not be working that well!! two nucleotides of TSO go to the UMI
    tsos <- c('ATAGACGAGATTGAAACTGCGCCCACCTATTAGCCACAGCCTATGCGTAGTAGGTATGTG',
              'CAAAGGCACATTGGGAGTCTAACCACAAGGGTCAGCGGTTCTATGCGTAGTAGGTATGTG',
              'CAAAGGCACAATGATGAGTTACCCACTCTCACGAAAACTCGTATGCGTAGTAGGTATGTT',
              'ATACTTAGGAAAGAACTCATTGCCCCTCCTCAATATGCCCTTATGCGTAGTAGGTATGTT',
              'GTTGTACCTTAGTGACACTCGAGAGACAAATGTATCGGAAATTTATTTATTTGTTTTTTT',
              'TNAGCGTTAAATGGCTAACTCACCACTGCATAGTACGTACTTTTATGCATAGTACCGCCA',
              'TNAACCCTGACCGTGATGTTCTCCAGACACTAGATAGACCGCACCTTTTTTTTTTTTTCT')

    ## fuzzy matching
    m <- aregexec(regex, c(test1, tsos), max.distance = list(substitutions = 4,
                                                             insertions = 0,
                                                             deletions = 0))
    foo <- do.call(rbind.data.frame, regmatches(c(test1, tsos), m))[,-1]
    colnames(foo) <- c('prepend', 'seq1', 'link1', 'seq2', 'link2', 'seq3', 'umi', 'tail')
    print(foo)


}


## ## extracts components of a fastq stanza and removes the prepended nts
## ## returns a tokenized fastq stanza, looking like this
## ##   prepend      seq1 link1      seq2 link2      seq3      umi              tail
## ## 1     TNA ATGTAATGG  GTGA ACGACCACC  GACA AAGGGAACT TCTTGATG    TTTTTTTTTTTTTT
## tokenize_stanza <- function(original_stanza, regex) {
##     x <- strcapture(pattern = regex,
##             x = original_stanza[2],
##             perl = TRUE,
##             proto = list(prepend = character(),
##                          seq1 = character(),
##                          link1 = character(),
##                          seq2 = character(),
##                          link2 = character(),
##                          seq3 = character(),
##                          umi = character(),
##                          tail = character()))

##     x$seq_trimmed <- paste(c(x$seq1, x$link1, x$seq2, x$link2, x$seq3, x$umi, x$tail), collapse = '')
##     x$qual_trimmed <- substr(original_stanza[4], nchar(x$prepend) + 1, nchar(original_stanza[4]))

##     return(x)
## }


## original_stanza <- c('@dT_prepend_2_nt_1_sub',
##                      'GCTGTACCTTAGTGAAATCCTGAAGACACAAGTTTCCATGCATGCTTTTTTTTTTTTTTTTT',
##                      '+',
##                      '12345678901234567890123456789012345678901234567890123456789012')

## original_stanza <- c('@dT_correct_prepend_1_nt_N_second_position',
##                      'ANGTACCTTAGTGAAATCCTGAAGACACAAGTTTCCATGCATGCTTTTTTTTTTTTTTTTTT',
##                      '+',
##                      '12345678901234567890123456789012345678901234567890123456789012')

## extracts components of a fastq stanza and removes the prepended nts
## returns a tokenized fastq stanza, looking like this
##   prepend      seq1 link1      seq2 link2      seq3      umi              tail
## 1     TNA ATGTAATGG  GTGA ACGACCACC  GACA AAGGGAACT TCTTGATG    TTTTTTTTTTTTTT
## mind we allow up to 4 of edit distance, now, as compared to the canonical regex
tokenize_stanza <- function(original_stanza, regex) {    
    ## fuzzy matching
    m <- aregexec(regex, original_stanza[2], max.distance = list(substitutions = 4,
                                                                 insertions = 0,
                                                                 deletions = 0))

    ## message(original_stanza[2])
    ## that is, doesn't match
    if (any(m[[1]] == '-1')) {
        x <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 8))
    } else { 
        x <- do.call(rbind.data.frame, regmatches(original_stanza[2], m))[,-1]
    }

    colnames(x) <- c('prepend', 'seq1', 'link1', 'seq2', 'link2', 'seq3', 'umi', 'tail')
    
    x$seq_trimmed <- paste(c(x$seq1, x$link1, x$seq2, x$link2, x$seq3, x$umi, x$tail), collapse = '')
    ## x$qual_trimmed <- substr(original_stanza[4], nchar(x$prepend) + 1, nchar(original_stanza[4]))
    x$qual_trimmed <- substr(original_stanza[4], start = nchar(original_stanza[4]) - nchar(x$seq_trimmed) + 1,
                             stop = nchar(original_stanza[4]))

    stopifnot(nchar(x$seq_trimmed) == nchar(x$qual_trimmed))

    return(x)
}


## original stanza is the original fastq, a string `\n` separated
## x is the result of capturing the regex groups, as generated by tokenize_stanza
## category the TSO/WTA/else nature of the read
reconstruct_stanza <- function(original_stanza, x, category = '') {
    trimmed_stanza <- sprintf('%s\n%s\n%s\n%s',
                              gsub('^@', paste0('@', category, '_'), original_stanza[1]),
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
has_Tn_in_tail <- function(x) {
      return(grepl(x = x$tail, pattern = 'TTTTT'))
}

## being x a tokenized result, as generated by tokenize_stanza
## currently (july 25th) its having the Tn subpattern in the tail,
##   as defined by has_Tn_in_tail()
is_Tn <- function(x) {
    return(has_Tn_in_tail(x))
}

## runs a normal grepl of character `pattern` on character `x`, and, if there is no match,
##   an agrepl (much slower) with `substitutions` substitutions
sequential_agrepl <- function(x, pattern, substitutions) {
    res <- grepl(x = x, pattern = pattern)
    if (!res & substitutions > 0) {
        res <- agrepl(x = x, pattern = pattern,
                  max.distance = list(substitutions = substitutions,
                                      insertions = 0,
                                      deletions = 0))
    }
    return(res)
}

## being x a tokenized result, as generated by tokenize_stanza
## returns a boolean: has a TSO or no
## it allows 1 subst
## this is under active development, the number of dels/subst/ins, or the TSO motif length,
##   are still being discussed
has_tso_in_tail <- function(x, substitutions = 1) {
    ## tso <- 'TATGCGTAGTAG'
    ## tso <- 'TATGCG'
    tso <- 'TATGCG'
    
    return(sequential_agrepl(x = x$tail, pattern = tso,
                  substitutions = substitutions))
}

## being x a tokenized result, as generated by tokenize_stanza
## currently (july 25th) its not having a prepended sequence; and having some TSO subpattern in the tail,
##   as defined by has_tso_in_tail()
is_tso <- function(x) {
    ## tryCatch({
    return(has_tso_in_tail(x) & nchar(x$prepend) == 0)
    ## }, error = function(y) print(dput(x)))
}
   
has_link_1 <- function(x, pattern, edit_dist = 1) {    
    return(sequential_agrepl(x = x$link1, pattern = pattern,
                  substitutions = edit_dist))
}

has_link_2 <- function(x, pattern, edit_dist = 1) {
    return(sequential_agrepl(x = x$link2, pattern = pattern,
                  substitutions = edit_dist))
}

## initialize_counts_object <- function() {
##     ## a log/summary to do some descriptive stats
##     counts <- list('non_matching' = 0,
##                    'GTGA' = list('GACA' = list('TSO' = 0,
##                                                'Tn' = as.list(setNames(rep(0, 20), sprintf('T%s', 1:20))),
##                                                'else' = 0)),
##                    'AATG' = list('CCAC' = list('TSO' = 0,
##                                                'Tn' = as.list(setNames(rep(0, 20), sprintf('T%s', 1:20))),
##                                                'else' = 0)),
##                    'other_linkers' = list('TSO' = 0,
##                                           'Tn' = as.list(setNames(rep(0, 20), sprintf('T%s', 1:20))),
##                                           'else' = 0))
##     return(counts)
## }

edit_cdna_id <- function(fastq_id, category){
    return(gsub('^@', paste0('@', category, '_'), fastq_id))
}

## original_stanza a barcode fastq stanza (4 lines, `\n` separeted)
## x is the tokenized result after capturing the regex groups
sort_reads_by_regex_match <- function(original_stanza, original_cdna, x) {

    
    ## Split reads into two buckets: TSO and not TSO.
    ## We do that according to the link1 and link2 sequences, and not checking whether the tail looks good or not
    ##   (even though TSO should have TSO seqs, and non TSO should have a polymer of Ts)

    ## the read does not match the regex, or doesn't have a tail
    if (any(is.na(x)) || x$tail == '') {
        ## message(sprintf('%s %s %s', x$linker1, x$linker2, x$tail))
        write(reconstruct_stanza(x = x, original_stanza = original_stanza, category = 'else'),
              file = stdout())

        write(edit_cdna_id(original_cdna, 'else'),
              file = stderr())
    }
    else {
        ## read is supposed to belong to bucket oligodT, but we check the tail
        if ((has_link_1(x = x, pattern = 'GTGA', edit_dist = 1) &
             has_link_2(x = x , pattern = 'GACA', edit_dist = 0)) |
            (has_link_1(x = x, pattern = 'GTGA', edit_dist = 0) &
             has_link_2(x = x , pattern = 'GACA', edit_dist = 1))) {
            
            if (is_Tn(x)) {
                write(reconstruct_stanza(x = x, original_stanza = original_stanza, category = 'WTA'),
                      file = stdout())

                write(edit_cdna_id(original_cdna, 'WTA'),
                      file = stderr())
                
            }
            else {
                write(reconstruct_stanza(x = x, original_stanza = original_stanza, category = 'else'),
                      file = stdout())

                write(edit_cdna_id(original_cdna, 'else'),
                      file = stderr())
            }
        }
        else if ((has_link_1(x = x, pattern = 'AATG', edit_dist = 1) &
                  has_link_2(x = x , pattern = 'CCAC', edit_dist = 0)) |
                 (has_link_1(x = x, pattern = 'AATG', edit_dist = 0) &
                  has_link_2(x = x , pattern = 'CCAC', edit_dist = 1))) {
            
            if (is_tso(x = x)) {                
                write(reconstruct_stanza(x = x, original_stanza = original_stanza, category = 'TSO'),
                      file = stdout())

                write(edit_cdna_id(original_cdna, 'TSO'),
                      file = stderr())

            }

            else {
                write(reconstruct_stanza(x = x, original_stanza = original_stanza, category = 'else'),
                      file = stdout())

                write(edit_cdna_id(original_cdna, 'else'),
                      file = stderr())
            }
        }
        else{
            ## the regex matches, but with noncanonical linkers; split by TSO/dT tails
            write(reconstruct_stanza(x = x, original_stanza = original_stanza, category = 'else'),
                  file = stdout())

             write(edit_cdna_id(original_cdna, 'else'),
                   file = stderr())
        }
    }
    return('done')
}


## iterates over each stanza from a fastq gz connection
process_stanza <- function(barcode_fn, cdna_fn, output_dir, regex) {
    ## the input filehandle
    fh <- gzfile(barcode_fn, 'r')
    cdna_fh <- gzfile(cdna_fn, 'r')
    

    
    while (TRUE) {
        original_barcode <- readLines(fh, n = 4)
        original_cdna <- readLines(cdna_fh, n = 4)

        
        if (length(original_barcode) == 0 | length(original_cdna) == 0) {
            break
        }

        ## tryCatch({
        counts <- sort_reads_by_regex_match(
            original_stanza = original_barcode,
            original_cdna = original_cdna,
            x = tokenize_stanza(original_stanza = original_barcode, regex = regex))
        ## }, error = function(x) message(original_barcode))
    }
    
    ## on.exit(close(fh))
    ## on.exit(close(cdna_fh))
}

process_stanza(barcode_fn = args$read1,
               cdna_fn = args$read2,
               regex = regex)

## cat(sprintf('End %s\t%s\n', args$output, Sys.time()))
