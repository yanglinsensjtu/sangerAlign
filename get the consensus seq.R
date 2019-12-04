source('select xstring by row name.R')
getAlignSeq <- function(seqset = xstringset,
                        geneid = geneid){
  
  for (i in seq_len(length(seqset)-1)) {
    i <- 1
    seqnames <- names(seqset)
    seqid <- seqnames[i]
    geneid <- seqnames[length(seqset)]
    
  
  a <- selectrow_byname(seqset,  
                        seqid,
                        geneid)
  alna <- msa(a,
              method = 'ClustalW',
              verbose = T, 
              order = 'input',
              gapOpening = 15,
              gapExtension = 6.66)
  consensus <- msaConsensusSequence(alna)
  startI <- str_locate_all(consensus,'[ATCG]')
  locate <- range(unlist(startI))
  tempseq <- toString(selectrow_byname(seqset,  
                                       seqid))
  seqset[seqid] <- str_sub(tempseq, start = locate[1], end = locate[2])
  }
  return(seqset)
}


