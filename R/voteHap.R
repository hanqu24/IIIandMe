#' Votting the final haplotype
#'
#' @param hap all possible haplotypes
#' @param npol pollen numbers
#' @param name names of pollens
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
#' VoteCount(hap, 5, colnames(input))

VoteCount <- function(hap, npol, name) {
  library(gtools)
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
    }
  }

  rownames(com) <- com$names
  com <- com[, -which(colnames(com) == "names")]
  dup <- com[!duplicated(lapply(com, summary))]

  count <- NULL
  for (i in 1:ncol(dup)) {
    nhap <- colSums(com == com[, match(colnames(dup)[i], colnames(com))])
    count <- c(count, sum(nhap == nrow(com)))
  }
  print("ratio is-----")
  print(count)
  print("----------")

  v1 <- dup[, which.max(count)]
  v2 <- flipFun(v1)
  vote <- cbind(v1, v2)
  rownames(vote) <- rownames(dup)
  colnames(vote) <- c("hap1", "hap2")

  if (ncol(dup) > 1) {
    tmp <- dup[, -which.max(count), drop = F]

    for (c in 1:ncol(tmp)) {
      whap <- colSums(com == com[, match(colnames(tmp)[c], colnames(com))])

      print("not major count is-----")
      se <- choice[which(whap == nrow(com)), ]
      if (!is.null(nrow(se))) {
        for (d in 1:nrow(se)) {
          print(name[se[d, ]])
        }
      } else {
        print(name[se])
      }

      print("----------")
    }
  }
  return(vote)
}
