library(sangerseqR)
library(Biostrings)
library(stringr)
library(msa)
# read in files -----------------------------------------------------------
sanger.resul.tpath <- '../sanger seq results/'
geneset <- readDNAStringSet(filepath = '../select genome for primer design.txt')
sgRNA <- readDNAStringSet('../sgRNAsequence.txt')
primer <- readDNAStringSet('../Primers_2019_11_19.fas')

# abi file rename ---------------------------------------------------------

library(stringr)

filelist <- dir(sanger.resul.tpath) %>% 
  str_subset(pattern = '\\.*.ab1$') %>% 
  str_sort(numeric = T)



for (i in seq_len(length(filelist))) {
  rename <- str_extract(filelist[i], '\\(.*') %>% 
    str_remove_all('\\)|\\[|\\]|\\(')
  filepath <- str_c(sanger.resul.tpath,filelist[i])
  filepathrename <- str_c(sanger.resul.tpath, rename)
  file.rename(from = filepath,to = filepathrename)
}


source('Non.specific.filter.R')#filter the nonspecific sequencing results
source('false discovery abi file filter.R')#filter the false discovery abi files
source('select xstring by row name.R')#select row of xstring object by name
source('get the consensus seq.R')
# read the abi files into the R env ---------------------------------------

filelist <- dir(sanger.resul.tpath) %>% 
  str_subset(pattern = '\\.*.ab1$') %>% str_sort(numeric = T)

genenames <- names(geneset)

for (i in seq_len(length(genenames))) {
  geneid <- genenames[i]


# remove false discovery abi files ----------------------------------------

FDabi.filter(geneid, 
             geneset,
             path = sanger.resul.tpath,
             consensus = 350)#remove the false discovery sequencing abi files

# after remove NF FD  -----------------------------------------------------


  filelist <- dir(sanger.resul.tpath) %>% 
    str_subset(pattern = '\\.*.ab1$') %>% 
    str_sort(numeric = T)
# tidy and import the primer into the R env -------------------------------



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
#seqset <- getAlignSeq(seqset)
#sgRNA <- DNAStringSet('GGAGGCCTGGGTTGCAGCGCTCCTGCAGATGGTAGAGGCGGCAGGGGCAGCTGAAGATGATCCCCTTCGGCGGACGGGGCGG')
#names(sgRNA) <- 'CD52sgRNA.b.site'
seqset <- append(seqset, sgRNA, after = length(seqset))
#seqset <- append(seqset, geneid.primer, after = length(seqset))


aln <- msa(seqset,
           method = 'ClustalOmega',
           verbose = T, 
           order = 'input')


# map the plot ------------------------------------------------------------
alnsgRNAstring <- toString(selectrow_byname(DNAStringSet(aln), names(sgRNA)))
locate <- str_locate(alnsgRNAstring, pattern = '[ATCG]+')
space <- round((100 - (locate[2]-locate[1]))/2)
ylim.a <- locate[1] - space
ylim.b <- locate[2] + space

print(aln, show="complete")
Biostrings::writeXStringSet(seqset,filepath = '../seqset.fasta')
source('msaprint.R')
msaprintPDF(aln,geneid, y = c(ylim.a, ylim.b))

}

source('change pdf file to png.R')
