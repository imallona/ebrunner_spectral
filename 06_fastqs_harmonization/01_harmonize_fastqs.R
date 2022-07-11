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
## $ Rscript 01_harmonize_fastqs.R -f test.fastq.gz -o gmoro
## [1] "Processing test.fastq.gz at gmoro"

## $ Rscript 01_harmonize_fastqs.R --help
## usage:01_harmonize_fastqs.R
##        [-h] [-f fastq] [-o fastq]

## optional arguments:
##   -h, --help            show this help message and exit
##   -f fastq, --fastq fastq
##                         Path to the R1 fastq file (gz compressed)
##   -o fastq, --output fastq
##                         Output folder name (creates the folder and overwrites content)

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

.libPaths('/home/imallona/R/R4_bioc314')

suppressPackageStartupMessages({
    library("argparse")
    library('R.utils')
})

parser <- ArgumentParser()

parser$add_argument("-f", "--fastq", type="character", default='none',
                    help="Path to the R1 fastq file (gz compressed)",
                    metavar="fastq")

parser$add_argument("-o", "--output", type="character", default='output',
                    help="Output folder name (creates the folder and overwrites content)",
                    metavar="fastq")

tryCatch({
    args <- parser$parse_args()
}, error = function(x) parser$print_help())


if ((args$fastq != 'none')) {
    print(sprintf('Processing %s at %s', args$fastq, args$output))
} else {
    parser$print_help()
}

dir.create(args$output, showWarnings = FALSE)

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

    strcapture(pattern = regex,
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
## returns a trimmed fastq stanza
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
## x is the result of capturing the regex groups
reconstruct_stanza <- function(original_stanza, x) {
    trimmed_stanza <- sprintf('%s\n%s\n%s\n%s\n',
                              original_stanza[1],
                              x$seq_trimmed,
                              original_stanza[3],
                              x$qual_trimmed)
    return(trimmed_stanza)
}

## being orignal_stanza a fastq stanza (4 lines, `\n` separeted)
## x is the result of capturing the regex groups
sort_reads_by_regex_match <- function(original_stanza, x, output_dir) {
    ## Split reads into two buckets: TSO and not TSO.
    ## We do that according to the link1 and link2 sequences, and not checking whether the tail looks good or not
    ##   (even though TSO should have TSO seqs, and non TSO should have a polymer of Ts)

    ## the read does not match the regex
    if (any(is.na(x))) {
        write(original_stanza, file = file.path(output_dir, 'non_matching.fastq'), append = TRUE)
    }
    else {
        if (x$link1 == 'GTGA' & x$link2 == 'GACA') {
            ## read belongs to bucket oligodT (we don't check the tail)
            write(reconstruct_stanza(x = x, original_stanza = original_stanza),
                  file = file.path(output_dir, 'GTGA_GACA_linkers.fastq'), append = TRUE)
        }
        else if (x$link1 == 'AATG' & x$link2 == 'CCAC') {
            ## read belongs to bucket TSO (we don't check the tail)
            write(reconstruct_stanza(x = x, original_stanza = original_stanza),
                  file = file.path(output_dir, 'AATG_CCAC_linkers.fastq'), append = TRUE)
        }
        else{
            ## read matches the regex but the two linkers are not as expected (?)
            ## this use case might not be needed
            write(original_stanza, file = file.path(output_dir, 'non_matching.fastq'), append = TRUE)
        }
    }
}

## iterates over each stanza from a fastq gz connection
process_stanzas <- function(fastq_gz, output_dir, regex) {
    ## the input filehandle
    fh <- gzfile(fastq_gz, 'r')

    ## the output 
    nonmatching <- file(file.path(output_dir, 'non_matching.fastq'), 'w')
    dt <- file(file.path(output_dir, 'GTGA_GACA_linkers.fastq'), 'w')
    tso <- file(file.path(output_dir, 'AATG_CCAC_linkers.fastq'), 'w')

    while (TRUE) {
        original_stanza <- readLines(fh, n = 4)
        if (length(original_stanza) == 0) {
            break
        }
        
        sort_reads_by_regex_match(original_stanza = original_stanza,
                                  x = tokenize_stanza(original_stanza = original_stanza, regex = regex),
                    output_dir = output_dir)
    }
    close(fh)

    ## compress the fastq outputs
    for (item in c('nonmatching', 'dt', 'tso')) {
        close(get(item))
    }
    
    ## compress the fastq outputs
    for (item in c('non_matching.fastq', 'GTGA_GACA_linkers.fastq', 'AATG_CCAC_linkers.fastq')) {
        ## close(file.path(output_dir, item))
        gzip(file.path(output_dir, item), ext = 'gz', overwrite = TRUE)
    }
}

process_stanzas(fastq_gz = args$fastq, output_dir = args$output, regex = regex)
