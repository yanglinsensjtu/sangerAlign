library(sangerseqR)
library(Biostrings)
library(stringr)
library(msa)
source('Non.specific.filter.R')#filter the nonspecific sequencing results
source('false discovery abi file filter.R')#filter the false discovery abi files
source('select xstring by row name.R')#select row of xstring object by name
source('get the consensus seq.R')

# read the abi files into the R env ---------------------------------------

sanger.resul.tpath <- '../sanger seq results/'
filelist <- dir(sanger.resul.tpath) %>% 
  str_subset(pattern = '\\.*.ab1$') %>% str_sort(numeric = T)

# read the sgRNA sequence into the R env ----------------------------------


sgRNA <- readDNAStringSet('../sgRNA sequence.fasta')



geneset <- readDNAStringSet(filepath = '../predict off target genes sequences.fasta')
genenames <- names(geneset)
geneid <- '8930'


FDabi.filter(geneid, geneset)#remove the false discovery sequencing abi files

# tidy and import the primer into the R env -------------------------------



primer <- readDNAStringSet('../Primers_2019_10_10.fas')
names(primer) <- str_extract(names(primer), pattern = '.*Primer..')
primernames <- str_subset(names(primer), pattern = geneid)
primernames <- str_sort(primernames)
geneid.primer <- selectrow_byname(primer, primernames)
geneid.primer[2] <- reverseComplement(geneid.primer[2])

# extract the sanger sequence into the DNSstringset obj -------------------


pat <- str_c('\\.*',geneid,'\\.*')
file <- str_subset(filelist, pattern = pat)
file <- str_sort(file, numeric = T)

seqset <- DNAStringSet(use.names=TRUE)
for (i in seq_len(length(file))) {
  filepath <- str_c(sanger.resul.tpath,file[i])
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
  seqset[i] <- reverseComplement(seqset[i])
  print(names(seqset[i]))
  print(sum1)
  print(sum2)
  if (sum1 > sum2) {
    next()
  }else{
    seqset[i] <- reverseComplement(seqset[i])
    if (str_detect(names(seqset[i]), 'F')) {
      names(seqset)[i] <- str_replace(names(seqset[i]), 'F', 'R')
    }else{
      names(seqset)[i] <- str_replace(names(seqset[i]), 'R', 'F')
    }
  }
}
seqset <- getAlignSeq(seqset)
seqset <- append(seqset, sgRNA, after = length(seqset))
seqset <- append(seqset, geneid.primer, after = length(seqset))


aln <- msa(seqset,
           method = 'ClustalW',
           verbose = T, 
           order = 'input')

print(aln, show="complete")

Biostrings::writeXStringSet(seqset,filepath = '../seqset.fasta')
source('msaprint.R')
msaprintPDF(aln,geneid)

# align the sgRNA and off target sequence ---------------------------------


offsetsgRNA <- DNAStringSet()
offsetsgRNA[1] <- toString(sgRNA)
offsetsgRNA[2] <- toString(geneset[3])
alnp <- msa(offsetsgRNA)
print(alnp, show = 'complete')



