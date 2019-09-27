library(sangerseqR)
library(Biostrings)
library(stringr)
library(muscle)

filelist <- dir('../../sanger 测序结果/') %>% 
  str_subset(pattern = '\\.*.ab1$')

seqset <- DNAStringSet()

geneid <- '3833'

pat <- str_c('\\.*',geneid,'\\.*')
file <- str_subset(filelist, pattern = '\\.*3833\\.*')
for (i in seq_len(length(file))) {
  filepath <- str_c('../../sanger 测序结果/',file[i])
  print(file[i])
  temp <- readsangerseq(filepath)
  seq <- as.character(temp@primarySeq)
  seqset[i] <- seq
}

aln <- muscle(seqset)
detail(aln)
