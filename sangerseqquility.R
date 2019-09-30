library(sangerseqR)
g50940.13 <- readsangerseq('../../sanger result/GENEID50940-13-M13F.ab1')
g50940.14 <- readsangerseq('../../sanger result/GENEID50940-14-M13F.ab1')
getTrace<- function(traceMatrix){
  traceMatrix <- tibble::as_tibble(traceMatrix(traceMatrix))
  colnames(traceMatrix) <- c('A','C','G','T')
  traceMatrix
}
g50940.13 <- getTrace(g50940.13)


g50940.14 <- getTrace(g50940.14)
score.sanger <- function(data){
 data <- mutate(data, score = (data$A<5)+
                               (data$C<5)+
                               (data$G <5)+ 
                               (data$`T`<5))
 data
}
g50940.13 <- score.sanger(g50940.13)
g50940.14 <- score.sanger(g50940.14)
