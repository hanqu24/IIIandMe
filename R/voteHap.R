#' Votting the final haplotype
#'
#' @param input Pollen data after pre-processing
#' @param hap_co haplotype & co lists
#' @param npol pollen numbers
#'
#' @return haplotype lists
#' @export
#'
#' @examples
#' library(IIIandme)
#' input<- PreProcessing(sample)
#' res<- HapCo(input, 5)
#' hap<- res[[1]]
#' co<- res[[2]]
#' VoteCount(sample, input, res, 5)

VoteCount <- function(sample, input, hap_co, npol) {
  library(gtools)
  hap<- hap_co[[1]]
  co<- hap_co[[2]]
  choice <- gtools::combinations(npol, 3, 5:(npol + 4))
  for (i in 1:length(hap)) {
    w <- which.max(c(sum(hap[[i]][, 1] == 0), sum(hap[[i]][, 2] == 0)))
    if (i == 1) {
      com <- data.frame(hap[[1]][, w])
      colnames(com) <- "hap_1"
      com$names <- rownames(hap[[1]])
    } else {
      com <- merge(com, data.frame(hap[[i]][, w], names = rownames(hap[[i]])), by = "names")
      colnames(com)[i + 1] <- paste("hap_", i, sep = "")
      com$names<- as.numeric(com$names)
      com<- com[order(com$names),]
    }
  }

  rownames(com) <- com$names
  com <- com[, -which(colnames(com) == "names")]
  dup <- com[!duplicated(lapply(com, summary))]

  count <- NULL
  whap<- list()
  for (i in 1:ncol(dup)) {
    nhap <- colSums(com == com[, match(colnames(dup)[i], colnames(com))])
    count <- c(count, sum(nhap == nrow(com)))
    whap[[i]]<- which(nhap == nrow(com))
  }
  print("ratio is-----")
  print(count)
  print("----------")

  v1 <- dup[, which.max(count)]
  v2 <- flipFun(v1)
  vote_hap <- cbind(v1, v2)
  rownames(vote_hap) <- rownames(dup)
  colnames(vote_hap) <- c("hap1", "hap2")

  wco<- whap[[which.max(count)]]
  aco<- do.call("rbind", co[c(wco)])
  nodup<- aggregate(list(numdup=rep(1,nrow(aco))), aco, length)
  nodup<- nodup[order(nodup$pollen, nodup$locus),]
  freq<- data.frame(table(choice[c(wco),]-4))
  freq$Var1<- paste('pol', freq$Var1, sep = '')
  nodup$percent<- nodup$numdup/freq[match(nodup$pollen, freq$Var1),2]
  nodup$vote<- rep('No', nrow(nodup))

  pol<- input[5:9]
  npol<- abs(pol[match(rownames(vote_hap), rownames(pol)),]-v1)

  for (a in 1: ncol(pol)) {
    pl <- crossoverCountFun(npol[, a])
    xx<- match(pl$locus, nodup$locus)
    if (prod(is.na(xx))!=1) {
      nodup[na.omit(xx), 7]<- 'Yes'
    }
  }

  #vote_co<- co[[as.numeric(strsplit(colnames(dup)[which.max(count)], '_')[[1]][2])]]
  vote_h<- apply(vote_hap, 2, function(x) ifelse(x == 1, sample$ref, sample$alt))
  vote<- list(vote_h, nodup)

  if (ncol(dup) > 1) {
    tmp <- dup[, -which.max(count), drop = F]

    for (c in 1:ncol(tmp)) {
      whap <- colSums(com == com[, match(colnames(tmp)[c], colnames(com))])

      print("not major count is-----")
      se <- choice[which(whap == nrow(com)), ]
      if (!is.null(nrow(se))) {
        for (d in 1:nrow(se)) {
          print(colnames(input)[se[d, ]])
        }
      } else {
        print(colnames(input)[se])
      }

      print("----------")
    }
  }
  return(vote)
}




