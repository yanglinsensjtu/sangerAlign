library(msa)
mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)
mySequences
myFirstAlignment <- msa(mySequences)
print(myFirstAlignment, show="complete")

bdir = getwd()
outdir = file.path(bdir,"alignements")

fasta_file = paste(outdir,".fasta",sep="")
output = file.path(outdir,"aln.tex", sep="")

dir.create(outdir,showWarnings = FALSE)
msaPrettyPrint(myFirstAlignment,
               output="tex", 
               showNames="none",
               showLogo="none", 
               askForOverwrite=FALSE, 
               verbose=FALSE,
               file = output, 
               alFile = fasta_file)
tools::texi2pdf(output, clean=TRUE)

