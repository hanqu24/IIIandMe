#' Plotting crossover positions
#'
#' @param sample Pollen raw data
#' @param chr chromosome name
#' @param co crossover table
#' @param scale data locus or real position
#' @param labels annotate label or not
#' @param label_font label font
#' @param label_angle label position
#'
#' @return crossover plot
#' @export
#'
#' @examples
#' library(IIIandMe)
#' input<- PreProcessing(sample)
#' res<- HapCo(input, 5)
#' vote<- VoteCount(sample, input, res, 5)
#' PlotCo(sample, chr='chr9', co=vote[[2]])

PlotCo<- function(sample, chr='chr9', co, scale='pos', export.options=F, labels=T,
                  label_font = 12,
                  label_angle = 0){
  library(chromoMap)
  chr<- chr
  anno<- data.frame()

  if (scale=='pos') {
    n<- sample[nrow(sample), 2]
    c<- data.frame(chr=colnames(sample)[6:ncol(sample)], start=c(1,1,1,1,1), end=n)
    co$pos<- as.numeric(co$pos)
    for (s in 1:length(unique(co$pollen))) {
      l<- which(co$pollen==unique(co$pollen)[s])
      if (length(l) %% 2 == 0) {
        x<- data.frame(matrix(co$pos[l], ncol = 2))
      }else{
        x<- data.frame(matrix(c(co$pos[l],n), ncol = 2))
      }
      colnames(x)<- c('start', 'end')
      name<- do.call(paste,c(chr,':', x[c('start')], '->', x[c('end')], sep=""))
      tmp<- data.frame(name= name, chr= rep(unique(co$pollen)[s], nrow(x)), x)
      anno<- rbind(anno, tmp)
    }
  }else{
    n<- nrow(sample)
    c<- data.frame(chr=colnames(sample)[6:ncol(sample)], start=c(1,1,1,1,1), end=n)
    co$locus<- as.numeric(co$locus)

    for (s in 1:length(unique(co$pollen))) {
      l<- which(co$pollen==unique(co$pollen)[s])
      if (length(l) %% 2 == 0) {
        x<- data.frame(matrix(co$locus[l], ncol = 2))
      }else{
        x<- data.frame(matrix(c(co$locus[l],n), ncol = 2))
      }
      colnames(x)<- c('start', 'end')
      name<- do.call(paste,c(chr,':', x[c('start')], '->', x[c('end')], sep=""))
      tmp<- data.frame(name= name, chr= rep(unique(co$pollen)[s], nrow(x)), x)
      anno<- rbind(anno, tmp)
    }
  }



  chromoMap(list(c),list(anno),segment_annotation = T, export.options=export.options,
            labels= labels,
            label_font = label_font,
            label_angle = label_angle)

}
