library(stringr)
sanger.resul.tpath <- '../../sanger result/'
filelist <- dir(sanger.resul.tpath) %>% 
  str_subset(pattern = '\\.*.ab1$')
geneid <- '19257'
pat <- str_c('\\.*',geneid,'\\.*')
file <- str_subset(filelist, pattern = pat)
file <- str_sort(file, numeric = T)
file
for (i in seq_len(length(file))) {
  filepath <- str_c(sanger.resul.tpath,file[i])
  rename <- str_replace(file[i], '19257', '197257')
  filepathrename <- str_c(sanger.resul.tpath, rename)
  file.rename(from = filepath,to = filepathrename)
}

