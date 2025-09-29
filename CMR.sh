Rscript - $1 $2 <<EOR
args=commandArgs(trailingOnly = TRUE)
bbj1=args[1];bbj2=args[2]
source('./CMR.r')
CMR(bbj1,bbj2,saveload='../MRresults')
EOR
