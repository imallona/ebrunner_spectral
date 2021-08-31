#!/usr/bin/env R
##
## Browses Spectral probes
##
## Started 31st Aug 2021
##
## Izaskun Mallona
## GPL v3


## Untested, based on https://bioconductor.org/packages/release/bioc/vignettes/epivizrChart/inst/doc/IntegrationWithShiny.html

library(epivizrChart)
library(rtracklayer)

epivizNav <- epivizNav(chr="chr11", start=118000000, end=121000000, interactive=TRUE)

## use local GTF?
genes_track <- epivizNav$add_genome(Homo.sapiens)

## probes?
probes <- rtracklayer::BEDFile("https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/refGene.hg19.bed.gz")

epiviz_igv <- epivizNav$plot(
                probes,
                datasource_name = "probes")


app <- shinyApp(
  ui=fluidPage(
    textInput('gene_loc', 'Enter Genomic Location (example: chr11:119000000 - 120000000', "chr11:118000000-121000000"),
    uiOutput("epivizChart")
  ),
  server=function(input, output, session) {
    
    renderEpiviz <- function() {
      output$epivizChart <- renderUI({
        epivizNav$render_component(shiny=TRUE)
      })
    }
    
    observeEvent(input$gene_loc, {
      loc <- input$gene_loc
      if(loc != "") {
        chr_split <- strsplit(loc, ":")
        chr <- chr_split[[1]][1]
        range_split <- strsplit(chr_split[[1]][2], "-")
        
        epivizNav$navigate(chr = chr, 
                           start = strtoi(range_split[[1]][1]), 
                           end = strtoi(range_split[[1]][2]))
      }
      renderEpiviz()
    })
    
    epivizNav$register_shiny_handler(session)
  }
)

app
