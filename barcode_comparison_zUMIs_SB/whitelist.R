### clearing environment 

rm(list=ls())

### loading barcodes

CLS1 <- read.table("~/RStudio/whitelist/CLS1", quote="\"", comment.char="")
CLS2 <- read.table("~/RStudio/whitelist/CLS2", quote="\"", comment.char="")
CLS3 <- read.table("~/RStudio/whitelist/CLS3", quote="\"", comment.char="")

CLS1 <- CLS1$V1
CLS2 <- CLS2$V1
CLS3 <- CLS3$V1

whitelist <- do.call(paste, expand.grid(CLS1, CLS2, CLS3, sep='', stringsAsFactors=FALSE))


write.table(whitelist, file = "~/RStudio/whitelist/whitelist.txt", sep = "\t",
            row.names = FALSE)

