library(sangerseqR)
library(Biostrings)
library(stringr)
library(msa)
filelist <- dir('../../sanger 测序结果/') %>% 
  str_subset(pattern = '\\.*.ab1$')

geneset <- readDNAStringSet(filepath = '../predict off target genes sequences.fasta')
geneset
geneid <- '144455'

pat <- str_c('\\.*',geneid,'\\.*')
file <- str_subset(filelist, pattern = pat)
file <- str_sort(file, numeric = T)

seqset <- DNAStringSet(use.names=TRUE)
for (i in seq_len(length(file))) {
  filepath <- str_c('../../sanger 测序结果/',file[i])
  print(file[i])
  temp <- readsangerseq(filepath)
  seq <- as.character(temp@primarySeq)
  seqset[i] <- seq
  names(seqset)[i] <- file[i]
}
seqset <- seqset[width(seqset) > 700]
seqset <- append(seqset,geneset[geneid], after = length(seqset))

for (i in seq_len(length(seqset)-1)) {
  msatemp<- msa(seqset[c(names(seqset)[i],geneid)],
                method = "ClustalOmega", 
                order = 'input')
  cn <- msaConsensusSequence(msatemp)
  sum1 <- sum(str_count(cn, pattern = 'A'),
              str_count(cn, pattern = 'T'),
              str_count(cn, pattern = 'C'),
              str_count(cn, pattern = 'G'))
  
  
  seqset[i] <- reverseComplement(seqset[i])
  msatemp<- msa(seqset[c(names(seqset)[i],geneid)],
                method = "ClustalOmega", 
                order = 'input')
  cn <- msaConsensusSequence(msatemp)
  sum2 <- sum(str_count(cn, pattern = 'A'),
              str_count(cn, pattern = 'T'),
              str_count(cn, pattern = 'C'),
              str_count(cn, pattern = 'G'))
  if (sum1 > sum2) {
    seqset[i] <- reverseComplement(seqset[i])
    next()
  }else{
    
    if (str_detect(names(seqset[i]), 'F')) {
      names(seqset)[i] <- str_replace(names(seqset[i]), 'F', 'R')
    }else{
      names(seqset)[i] <- str_replace(names(seqset[i]), 'R', 'F')
    }
   
  }
}
aln <- msa(seqset,
           method = "ClustalOmega", 
           verbose = T, 
           order = 'input')

print(aln, show="complete")

Biostrings::writeXStringSet(seqset,filepath = '../seqset.fasta')


