#!/usr/bin/env/R

## tso <- 'TATGCGTAGTAGGTATG'
## dt <- paste(rep('T', 25), collapse = '')

regex<- '([ACTGN]{0,3})([ACTGN]{9})(GTGA|AATG{1})([ACTGN]{9})(GACA|CCAC){1}([ACTGN]{9})([ACTGN]{8})(.*)'



test1 <- c('TNAATGTAATGGGTGAACGACCACCGACAAAGGGAACTTCTTGATGTTTTTTTTTTTTTT',
           'CNCATTGCAAATGAATCCTGAACCACTTCAGCTCAGTTAAATATATGCGTAGGAGGTATG',
           'GNACACACAAAGTGAAGATAGTTCGACAACACCTTAGCTGCTTACTTTTTTTGTTTTTTT',
           'ANTTAGATGAATGCAGAAATCGCCACGATGGTCCACAGTTTCTTATGCCTACAACCGATG')

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
