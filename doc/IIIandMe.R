## ----install, eval=FALSE, include=FALSE---------------------------------------
#  ### Install dependencies ahead
#  install.packages('devtools')
#  
#  devtools::install_github('Jialab-UCR/IIIandMe')

## ----raw, echo=TRUE-----------------------------------------------------------
### Load sample data

library(IIIandMe)
data(sample)
head(sample)


## ----pre, echo=TRUE-----------------------------------------------------------

input<- PreProcessing(sample)
head(input)


## ----res, include=FALSE-------------------------------------------------------

res<- HapCo(input, 5, filterGenoError = F)
hap<- res[[1]]
co<- res[[2]]

## ----hap, echo=TRUE-----------------------------------------------------------
head(hap[[1]])
co[[1]]

## ----vote, echo=TRUE----------------------------------------------------------
vote<- VoteCount(sample, input, res, 5)
hap<- vote[[1]]
co<- vote[[2]]

## ----show---------------------------------------------------------------------
head(hap)
co

## ----plot, echo=TRUE----------------------------------------------------------
PlotCo(sample, chr='chr9', co=vote[[2]])

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

