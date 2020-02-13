library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg18)
library(tidyverse)

bval_path <- "bkgd/"
chain_path <- "hg18ToHg19.over.chain"

# load all the original bkgd files from McVicker et al.
bval_files <- list.files(bval_path, full.names=TRUE, ".*.bkgd")
bval_df_list <- lapply(bval_files, function(filename) {
    read.table(filename, col.names=c("bval", "length")) %>%
        mutate(chr=str_replace(basename(filename), ".bkgd", ""),
               end=cumsum(length),
               start=c(1, (end + 1)[-n()])) %>%
        select(chr, start, end, bval)
})

## convert the list of dataframes into a GRanges hg18 object
bval_regions_hg18 <- bind_rows(bval_df_list) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
seqinfo(bval_regions_hg18) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)

## perform the liftOver from hg18 to hg19
chain <- import.chain(chain_path)
bval_regions <- liftOver(bval_regions_hg18, chain) %>% unlist

# save output to bed
as.data.frame(bval_regions) %>%
    mutate(start = start - 1) %>%
    select(chrom = seqnames, start, end, bval) %>%
    write_tsv("bvalues_hg19.bed.gz", col_names = T)
