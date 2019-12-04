#filter the good sequencing abi but false discovery file
library(sangerseqR)
library(Biostrings)
library(stringr)
library(msa)
source('sangerseqquility.R')
FDabi.filter <- function(geneid = geneid,
                         geneset = geneset,
                         path = sanger.resul.tpath,
                         consensus = numeric){
  filelist <- dir(sanger.resul.tpath) %>% 
    str_subset(pattern = '\\.*.ab1$') %>% 
    str_sort(numeric = T)
  FD <- file.path(path, 'FD')
  if (file.exists(FD)) {
    print('The FD folder already existed')
  }else{
    dir.create(FD)
  }
  
  filelist <- str_subset(filelist, pattern = geneid)
  for (i in seq_len(length(filelist))) {
    seqset <- DNAStringSet(use.names=TRUE)
    filepath <- str_c(path,filelist[i])
    if(file.exists(filepath)){
      temp <- readsangerseq(filepath)
    }else{
      print(str_c(filelist[i],'not existed'))
      next()
    }
    
    seq <- as.character(temp@primarySeq)
    seqset[1] <- seq
    names(seqset)[1] <- as.character(filelist[i])
    seqset <- append(seqset, geneset[geneid], after = length(seqset))
    seqset
    msatemp<- msa(seqset,
                  method = "ClustalW", 
                  order = 'input')
    cn1 <- msaConsensusSequence(msatemp)
    sum1 <- sum(str_count(cn1, pattern = 'A'),
                str_count(cn1, pattern = 'T'),
                str_count(cn1, pattern = 'C'),
                str_count(cn1, pattern = 'G'))
    seqset[1] <- reverseComplement(seqset[1])
    msatemp<- msa(seqset,
                  method = "ClustalW", 
                  order = 'input')
    cn2 <- msaConsensusSequence(msatemp)
    sum2 <- sum(str_count(cn2, pattern = 'A'),
                str_count(cn2, pattern = 'T'),
                str_count(cn2, pattern = 'C'),
                str_count(cn2, pattern = 'G'))
    if (sum1 < consensus && sum2 < consensus) {
      print(paste(filelist[i],'is false discovery abi file'))
      file.copy(filepath, FD)
      file.remove(filepath)
    }else{
      next()
    }
  }
  geneset
}
